#include <fstream>
#include <iostream>
#include <iomanip>
#include "src/model/stoch_oil/StochOilMethod.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using namespace stoch_oil;

StochOilMethod::StochOilMethod(Model* _model) : AbstractMethod<Model>(_model)
{
	const int strNum0 = var_size0 * model->cellsNum;
	y0 = new double[strNum0];
	ind_i0 = new int[Mesh::stencil * var_size0 * strNum0];
	ind_j0 = new int[Mesh::stencil * var_size0 * strNum0];
	//cols = new int[strNum];
	a0 = new double[Mesh::stencil * var_size0 * strNum0];
	ind_rhs0 = new int[strNum0];
	rhs0 = new double[strNum0];

	const int strNum1 = var_size1 * model->cellsNum;
	y1 = new double[strNum1];
	ind_i1 = new int[Mesh::stencil * var_size1 * strNum1];
	ind_j1 = new int[Mesh::stencil * var_size1 * strNum1];
	//cols = new int[strNum];
	a1 = new double[Mesh::stencil * var_size1 * strNum1];
	ind_rhs1 = new int[strNum1];
	rhs1 = new double[strNum1];

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
	delete[] y0, y1;
	//for (int i = 0; i < Model::var_size * size; i++)
	//	delete[] jac[i];
	//delete[] jac;

	delete[] ind_i0, ind_j0, ind_rhs0;
	delete[] a0, rhs0;

	delete[] ind_i1, ind_j1, ind_rhs1;
	delete[] a1, rhs1;

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

		for (size_t i = 0; i < var_size0; i++)
			for (const int idx : cell.stencil)
				for (size_t j = 0; j < var_size0; j++)
				{
					ind_i0[counter] = var_size0 * cell.id + i;			ind_j0[counter++] = var_size0 * idx + j;
				}
	}

	elemNum0 = counter;

	for (int i = 0; i < var_size0 * model->cellsNum; i++)
		ind_rhs0[i] = i;
};
void StochOilMethod::fillIndices1()
{
	int counter = 0;

	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = mesh->cells[i];
		getMatrixStencil(cell);

		for (size_t i = 0; i < var_size1; i++)
			for (const int idx : cell.stencil)
				for (size_t j = 0; j < var_size1; j++)
				{
					ind_i1[counter] = var_size1 * cell.id + i;			ind_j1[counter++] = var_size1 * idx + j;
				}
	}

	elemNum1 = counter;

	for (int i = 0; i < var_size1 * model->cellsNum; i++)
		ind_rhs1[i] = i;
};
void StochOilMethod::start()
{
	step_idx = 0;

	fillIndices0();
	solver0.Init(var_size0 * model->cellsNum, 1.e-15, 1.e-15);
	fillIndices1();
	solver1.Init(var_size1 * model->cellsNum, 1.e-15, 1.e-15);

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
		solver0.Assemble(ind_i0, ind_j0, a0, elemNum0, ind_rhs0, rhs0);
		solver0.Solve(PRECOND::ILU_SIMPLE);
		copySolution0(solver0.getSolution());

		err_newton = convergance(cellIdx, varIdx);
		averValue(averVal);
		for (int i = 0; i < var_size0; i++)
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
	{
		auto& var = (*model)[i].u_next0;
		var.p0 += sol[Model::var_size0 * i];
	}
}
void StochOilMethod::copySolution1(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < size; i++)
	{
		auto& var = (*model)[i].u_next1;
		var.p2 += sol[Model::var_size1 * i];
		var.Cfp += sol[Model::var_size1 * i + 1];
	}
}
void StochOilMethod::computeJac0()
{
	trace_on(0);

	for (size_t i = 0; i < size; i++)
		model->x0[i].p0 <<= model->u_next0[i * var_size0];

	for (int i = 0; i < size; i++)
	{
		const auto& cell = mesh->cells[i];

		if (cell.type == elem::QUAD)
			model->h0[i * var_size0] = model->solveInner0(cell).p0;
		else if (cell.type == elem::BORDER)
			model->h0[i * var_size0] = model->solveBorder0(cell).p0;
	}

	for (const auto& well : model->wells)
	{
		if (well.cur_bound)
			model->h0[well.cell_id * var_size0] += model->solveSource0(well).p0;
	}

	for (int i = 0; i < var_size0 * size; i++)
		model->h0[i] >>= y0[i];

	trace_off();
}
void StochOilMethod::computeJac1()
{
	trace_on(1);

	for (size_t i = 0; i < size; i++)
	{
		model->x1[i].p2 <<= model->u_next1[i * var_size1];
		model->x1[i].Cfp <<= model->u_next1[i * var_size1 + 1];
	}

	for (int i = 0; i < size; i++)
	{
		const auto& cell = mesh->cells[i];

		if (cell.type == elem::QUAD)
		{
			const auto& system = model->solveInner1(cell);
			model->h1[i * var_size1] = system.p2;
			model->h1[i * var_size1 + 1] = system.Cfp;
		}
		else if (cell.type == elem::BORDER)
		{
			const auto& system = model->solveBorder1(cell);
			model->h1[i * var_size1] = system.p2;
			model->h1[i * var_size1 + 1] = system.Cfp;
		}
	}

	for (const auto& well : model->wells)
	{
		if (well.cur_bound)
		{
			const auto& system = model->solveSource1(well);
			model->h1[well.cell_id * var_size1] += system.p2;
			model->h1[well.cell_id * var_size1 + 1] += system.Cfp;
		}
	}

	for (int i = 0; i < var_size1 * size; i++)
		model->h1[i] >>= y1[i];

	trace_off();
}
void StochOilMethod::fill0()
{
	sparse_jac(0, var_size0 * model->cellsNum, var_size0 * model->cellsNum, repeat,
		&model->u_next0[0], &elemNum0, (unsigned int**)(&ind_i0), (unsigned int**)(&ind_j0), &a0, options);

	int counter = 0;
	for (int j = 0; j < size; j++)
	{
		const auto& cell = mesh->cells[j];
		//getMatrixStencil(cell);
		for (int i = 0; i < var_size0; i++)
		{
			const int str_idx = var_size0 * cell.id + i;
			/*for (const int idx : stencil_idx)
			{
			for (int j = 0; j < Model::var_size; j++)
			a[counter++] = jac[str_idx][Model::var_size * idx + j];
			}*/

			rhs0[str_idx] = -y0[str_idx];
		}
		//stencil_idx.clear();
	}
}
void StochOilMethod::fill1()
{
	sparse_jac(1, var_size1 * model->cellsNum, var_size1 * model->cellsNum, repeat,
		&model->u_next1[0], &elemNum1, (unsigned int**)(&ind_i1), (unsigned int**)(&ind_j1), &a1, options);

	int counter = 0;
	for (int j = 0; j < size; j++)
	{
		const auto& cell = mesh->cells[j];
		//getMatrixStencil(cell);
		for (int i = 0; i < var_size1; i++)
		{
			const int str_idx = var_size1 * cell.id + i;
			/*for (const int idx : stencil_idx)
			{
			for (int j = 0; j < Model::var_size; j++)
			a[counter++] = jac[str_idx][Model::var_size * idx + j];
			}*/

			rhs1[str_idx] = -y1[str_idx];
		}
		//stencil_idx.clear();
	}
}