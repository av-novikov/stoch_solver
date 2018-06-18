#include <fstream>
#include <iostream>
#include <iomanip>
#include "src/model/oil/OilMethod.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using namespace oil;

OilMethod::OilMethod(Model* _model) : AbstractMethod<Model>(_model)
{
	const int strNum = var_size * model->cellsNum;

	y = new double[strNum];
	ind_i = new int[Mesh::stencil * var_size * strNum];
	ind_j = new int[Mesh::stencil * var_size * strNum];
	//cols = new int[strNum];
	a = new double[Mesh::stencil * var_size * strNum];
	ind_rhs = new int[strNum];
	rhs = new double[strNum];

	//options[0] = 0;          /* sparsity pattern by index domains (default) */
	//options[1] = 0;          /*                         safe mode (default) */
	//options[2] = 0;          /*              not required if options[0] = 0 */
	//options[3] = 0;          /*                column compression (default) */

	pvd.open("snaps/Oil.pvd", std::ofstream::out);
	pvd << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd << "\t<Collection>\n";

	plot_P.open("snaps/P.dat", std::ofstream::out);
	plot_Q.open("snaps/Q.dat", std::ofstream::out);
};
OilMethod::~OilMethod()
{
	delete[] y;
	//for (int i = 0; i < Model::var_size * size; i++)
	//	delete[] jac[i];
	//delete[] jac;

	delete[] ind_i, ind_j, ind_rhs;
	//delete[] cols;
	delete[] a, rhs;

	plot_P.close();
	plot_Q.close();
	pvd << "\t</Collection>\n";
	pvd << "</VTKFile>\n";
	pvd.close();
};
void OilMethod::writeData()
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
		"0\" file=\"Oil_" + to_string(step_idx) + ".vtu\"/>\n";
}
void OilMethod::control()
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
void OilMethod::start()
{
	step_idx = 0;

	fillIndices();
	solver.Init(Model::var_size * model->cellsNum, 1.e-15, 1.e-15);

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
void OilMethod::solveStep()
{
	int cellIdx, varIdx, iterations;
	double err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	
	iterations = 0;
	while (err_newton > 1.e-4 && dAverVal[0] > 1.e-7 && iterations < 20)
	{
		copyIterLayer();
		computeJac();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve(PRECOND::ILU_SIMPLE);
		copySolution(solver.getSolution());

		err_newton = convergance(cellIdx, varIdx);
		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;
		iterations++;
	}

	std::cout << "Newton Iterations = " << iterations << std::endl;
}
void OilMethod::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < size; i++)
	{
		auto& var = (*model)[i].u_next;
		var.p0 += sol[Model::var_size * i];
	}
}
void OilMethod::computeJac()
{
	trace_on(0);

	for (size_t i = 0; i < size; i++)
		model->x[i].p0 <<= model->u_next[i * var_size];

	for (int i = 0; i < size; i++)
	{
		const auto& cell = mesh->cells[i];

		if (cell.type == elem::QUAD)
			model->h[i] = model->solveInner(cell);
		else if (cell.type == elem::BORDER)
			model->h[i] = model->solveBorder(cell);
	}

	for (const auto& well : model->wells)
	{
		if (well.cur_bound)
			model->h[well.cell_id] += well.cur_rate * model->props_oil.rho_stc / mesh->hz * model->ht;
	}

	for (int i = 0; i < Model::var_size * size; i++)
		model->h[i] >>= y[i];

	trace_off();
}
void OilMethod::fill()
{
	sparse_jac(0, Model::var_size * model->cellsNum, Model::var_size * model->cellsNum, repeat,
		&model->u_next[0], &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	int counter = 0;
	for (int j = 0; j < size; j++)
	{
		const auto& cell = mesh->cells[j];
		//getMatrixStencil(cell);
		for (int i = 0; i < Model::var_size; i++)
		{
			const int str_idx = Model::var_size * cell.id + i;
			/*for (const int idx : stencil_idx)
			{
			for (int j = 0; j < Model::var_size; j++)
			a[counter++] = jac[str_idx][Model::var_size * idx + j];
			}*/

			rhs[str_idx] = -y[str_idx];
		}
		//stencil_idx.clear();
	}
}