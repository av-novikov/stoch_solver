#include <fstream>
#include <iostream>
#include <iomanip>
#include "src/model/stoch_oil/StochOilMethod.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using namespace stoch_oil;

StochOilMethod::StochOilMethod(Model* _model) : AbstractMethod<Model>(_model)
{
	const int strNum0 = model->cellsNum;
	y_p0 = new double[strNum0];
	ind_i_p0 = new int[Mesh::stencil * strNum0];
	ind_j_p0 = new int[Mesh::stencil * strNum0];
	//cols = new int[strNum];
	a_p0 = new double[Mesh::stencil * strNum0];
	ind_rhs_p0 = new int[strNum0];
	rhs_p0 = new double[strNum0];

	const int strNum1 = model->cellsNum;
	y_Cfp = new double[strNum1];
	ind_i_Cfp = new int[Mesh::stencil * strNum1];
	ind_j_Cfp = new int[Mesh::stencil * strNum1];
	//cols = new int[strNum];
	a_Cfp = new double[Mesh::stencil * strNum1];
	ind_rhs_Cfp = new int[strNum1];
	rhs_Cfp = new double[strNum1];

	//options[0] = 0;          /* sparsity pattern by index domains (default) */
	//options[1] = 0;          /*                         safe mode (default) */
	//options[2] = 0;          /*              not required if options[0] = 0 */
	//options[3] = 0;          /*                column compression (default) */

	pvd.open("snaps/StochOil.pvd", std::ofstream::out);
	pvd << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd << "\t<Collection>\n";

	plot_P.open("snaps/P.dat", std::ofstream::out);
	plot_Q.open("snaps/Q.dat", std::ofstream::out);
};
StochOilMethod::~StochOilMethod()
{
	delete[] y_p0, y_Cfp;
	//for (int i = 0; i < Model::var_size * size; i++)
	//	delete[] jac[i];
	//delete[] jac;

	delete[] ind_i_p0, ind_j_p0, ind_rhs_p0;
	delete[] a_p0, rhs_p0;

	delete[] ind_i_Cfp, ind_j_Cfp, ind_rhs_Cfp;
	delete[] a_Cfp, rhs_Cfp;

	plot_P.close();
	plot_Q.close();
	pvd << "\t</Collection>\n";
	pvd << "</VTKFile>\n";
	pvd.close();
};
void StochOilMethod::writeData()
{
	double p = 0.0, q = 0;

	plot_Q << cur_t * t_dim / 3600.0;
	plot_P << cur_t * t_dim / 3600.0;

	for (const auto& well : model->wells)
	{
		plot_Q << "\t" << model->getRate(well) * model->Q_dim * 86400.0;
		plot_P << "\t" << model->getPwf(well) * model->P_dim / BAR_TO_PA;
	}

	plot_Q << std::endl;
	plot_P << std::endl;

	pvd << "\t\t<DataSet part=\"0\" timestep=\"" + to_string(cur_t * t_dim / 3600.0) +
		"0\" file=\"StochOil_" + to_string(step_idx) + ".vtu\"/>\n";
}
void StochOilMethod::control()
{
	writeData();

	if (cur_t >= model->wells[0].period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	model->ht *= 1.5;
	//if (model->ht <= model->ht_max && iterations < 6)
	//	model->ht = model->ht * 1.5;
	//else if (iterations > 6 && model->ht > model->ht_min)
	//	model->ht = model->ht / 1.5;

	if (cur_t + model->ht > model->wells[0].period[curTimePeriod])
		model->ht = model->wells[0].period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void StochOilMethod::fillIndices0()
{
	int counter = 0;

	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = mesh->cells[i];
		getMatrixStencil(cell);

		for (size_t i = 0; i < var_size; i++)
			for (const int idx : cell.stencil)
				for (size_t j = 0; j < var_size; j++)
				{
					ind_i_p0[counter] = var_size * cell.id + i;			ind_j_p0[counter++] = var_size * idx + j;
				}
	}

	elemNum_p0 = counter;

	for (int i = 0; i < var_size * model->cellsNum; i++)
		ind_rhs_p0[i] = i;
};
void StochOilMethod::fillIndices1()
{
	int counter = 0;

	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = mesh->cells[i];
		getMatrixStencil(cell);

		for (size_t i = 0; i < var_size; i++)
			for (const int idx : cell.stencil)
				for (size_t j = 0; j < var_size; j++)
				{
					ind_i_Cfp[counter] = var_size * cell.id + i;			ind_j_Cfp[counter++] = var_size * idx + j;
				}
	}

	elemNum_Cfp = counter;

	for (int i = 0; i < var_size * model->cellsNum; i++)
		ind_rhs_Cfp[i] = i;
};
void StochOilMethod::start()
{
	step_idx = 0;

	fillIndices0();
	solver0.Init(model->cellsNum, 1.e-15, 1.e-15);
	fillIndices1();
	solver1.Init(model->cellsNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);
	while (cur_t < Tt)
	{
		control();
		model->snapshot_all(step_idx++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << std::endl;
		cout << std::setprecision(6);
		cout << "time = " << cur_t << std::endl;
	}

	model->snapshot_all(step_idx);
	writeData();
}
void StochOilMethod::solveStep()
{
	int cellIdx, varIdx, iterations;
	double err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	
	iterations = 0;
	while (err_newton > 1.e-4 && dAverVal[0] > 1.e-7 && iterations < 20)
	{
		copyIterLayer();
		computeJac0();
		fill0();
		solver0.Assemble(ind_i_p0, ind_j_p0, a_p0, elemNum_p0, ind_rhs_p0, rhs_p0);
		solver0.Solve(PRECOND::ILU_SIMPLE);
		copySolution0(solver0.getSolution());

		err_newton = convergance(cellIdx, varIdx);
		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		/*computeJac1();
		fill1();
		solver1.Assemble(ind_i1, ind_j1, a1, elemNum1, ind_rhs1, rhs1);
		solver1.Solve(PRECOND::ILU_SIMPLE);
		copySolution1(solver1.getSolution());*/

		iterations++;
	}

	std::cout << "Newton Iterations = " << iterations << std::endl;
}
void StochOilMethod::copySolution0(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < size; i++)
		model->p0_next[i] += sol[i];
}
void StochOilMethod::copySolution1(const int cell_id, const paralution::LocalVector<double>& sol)
{
	for (int i = cell_id * size; i < (cell_id + 1) * size; i++)
		model->Cfp_next[i] += sol[i];
}
void StochOilMethod::computeJac0()
{
	trace_on(0);

	for (size_t i = 0; i < size; i++)
		model->x_p0[i] <<= model->p0_next[i * var_size];

	for (int i = 0; i < size; i++)
	{
		const auto& cell = mesh->cells[i];

		if (cell.type == elem::QUAD)
			model->h_p0[i * var_size] = model->solveInner0(cell);
		else if (cell.type == elem::BORDER)
			model->h_p0[i * var_size] = model->solveBorder0(cell);
	}

	for (const auto& well : model->wells)
		if (well.cur_bound)
			model->h_p0[well.cell_id * var_size] += model->solveSource0(well);

	for (int i = 0; i < var_size * size; i++)
		model->h_p0[i] >>= y_p0[i];

	trace_off();
}
void StochOilMethod::computeJac1()
{
	trace_on(1);

	for (size_t i = 0; i < size * size; i++)
		model->x_Cfp[i] <<= model->Cfp_next[i * var_size];

	for (int i = 0; i < size * size; i++)
	{
		const auto& cell = mesh->cells[i];

		if (cell.type == elem::QUAD)
			model->h_Cfp[i * var_size] = model->solveInner1(cell);
		else if (cell.type == elem::BORDER)
			model->h_Cfp[i * var_size] = model->solveBorder1(cell);
	}

	for (const auto& well : model->wells)
		if (well.cur_bound)
			model->h_Cfp[well.cell_id * var_size] += model->solveSource1(well);

	for (int i = 0; i < size * size; i++)
		model->h_Cfp[i] >>= y_Cfp[i];

	trace_off();
}
void StochOilMethod::fill0()
{
	sparse_jac(0, model->cellsNum, model->cellsNum, repeat,
		&model->p0_next[0], &elemNum_p0, (unsigned int**)(&ind_i_p0), (unsigned int**)(&ind_j_p0), &a_p0, options);

	int counter = 0;
	for (int j = 0; j < size; j++)
	{
		const auto& cell = mesh->cells[j];
		rhs_p0[cell.id] = -y_p0[cell.id];
	}
}
void StochOilMethod::fill1()
{
	sparse_jac(1, model->cellsNum, model->cellsNum, repeat,
		&model->Cfp_next[0], &elemNum_Cfp, (unsigned int**)(&ind_i_Cfp), (unsigned int**)(&ind_j_Cfp), &a_Cfp, options);

	int counter = 0;
	for (int j = 0; j < size; j++)
	{
		const auto& cell = mesh->cells[j];
		rhs_Cfp[cell.id] = -y_Cfp[cell.id];
	}
}