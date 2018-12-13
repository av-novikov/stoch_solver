#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkQuad.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkHexahedron.h>

#include "src/utils/VTKSnapshotter.hpp"
#include "src/utils/utils.h"

#include "src/model/oil/Oil.hpp"
#include "src/model/stoch_oil/StochOil.hpp"

using namespace std;
using namespace snapshotter;

template<class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter(const Model* _model) : model(_model), mesh(_model->getMesh())
{
	R_dim = model->R_dim;
	pattern = prefix + "Mesh_%{STEP}.vtu";
}
VTKSnapshotter<oil::Oil>::VTKSnapshotter(const oil::Oil* _model) : model(_model), mesh(_model->getMesh())
{
	R_dim = model->R_dim;
	pattern = prefix + "Oil_%{STEP}.vtu";

	num_x = mesh->num_x;	num_y = mesh->num_y;
}
VTKSnapshotter<stoch_oil::StochOil>::VTKSnapshotter(const stoch_oil::StochOil* _model) : model(_model), mesh(_model->getMesh())
{
	R_dim = model->R_dim;
	pattern = prefix + "StochOil_%{STEP}.vtu";

	num_x = mesh->num_x;	num_y = mesh->num_y;
}
template<class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
}
template<class modelType>
string VTKSnapshotter<modelType>::replace(string filename, string from, string to)
{
	size_t start_pos = 0;
	while ((start_pos = filename.find(from, start_pos)) != string::npos)
	{
		filename.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return filename;
}
template<class modelType>
std::string VTKSnapshotter<modelType>::getFileName(const int i)
{
	string filename = pattern;
	return replace(filename, "%{STEP}", to_string(i));
}

template<class modelType>
void VTKSnapshotter<modelType>::dump(const int i)
{
}
void VTKSnapshotter<oil::Oil>::dump(const int snap_idx)
{
	using namespace oil;
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	auto p0 = vtkSmartPointer<vtkDoubleArray>::New();
	p0->SetName("p0");

	points->Allocate((num_x + 1) * (num_y + 1));
	cells->Allocate(num_x * num_y);
	double hx = mesh->hx / (double)num_x;
	double hy = mesh->hy / (double)num_y;

	for (int i = 0; i < num_x + 1; i++)
		for (int j = 0; j < num_y + 1; j++)
		{
			const Cell& cell = mesh->cells[(num_y + 2) * i + j];
			points->InsertNextPoint(R_dim * (cell.cent.x + cell.hx / 2), R_dim * (cell.cent.y + cell.hy / 2), 0.0);
		}
	grid->SetPoints(points);

	size_t x_ind, y_ind;
	for (const auto& cell : mesh->cells)
	{
		if (cell.type == elem::QUAD)
		{
			x_ind = cell.id / (num_y + 2) - 1;
			y_ind = cell.id % (num_y + 2) - 1;

			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
			quad->GetPointIds()->SetId(0, y_ind + x_ind * (num_y + 1));
			quad->GetPointIds()->SetId(1, y_ind + x_ind * (num_y + 1) + 1);
			quad->GetPointIds()->SetId(2, y_ind + (x_ind + 1) * (num_y + 1) + 1);
			quad->GetPointIds()->SetId(3, y_ind + (x_ind + 1) * (num_y + 1));
			cells->InsertNextCell(quad);

			const auto& var = (*model)[cell.id].u_next;
			p0->InsertNextValue(var.p0 * model->P_dim / BAR_TO_PA);
		}
	}
	grid->SetCells(VTK_QUAD, cells);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(p0);

	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<stoch_oil::StochOil>::dump(const int snap_idx)
{
	using namespace stoch_oil;
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto cells = vtkSmartPointer<vtkCellArray>::New();
	auto p0 = vtkSmartPointer<vtkDoubleArray>::New();
	p0->SetName("p0");
	auto p2 = vtkSmartPointer<vtkDoubleArray>::New();
	p2->SetName("p2");
	
    std::vector<vtkSmartPointer<vtkDoubleArray>> pwf_perm_corr;
    std::vector<vtkSmartPointer<vtkDoubleArray>> pwf_pres_corr;
    std::vector<vtkSmartPointer<vtkDoubleArray>> q_perm_corr;
    std::vector<vtkSmartPointer<vtkDoubleArray>> q_pres_corr;
    std::vector<vtkSmartPointer<vtkDoubleArray>> Cf_well;
    pwf_perm_corr.resize(model->wells.size());
    pwf_pres_corr.resize(model->wells.size());
    q_perm_corr.resize(model->wells.size());
    q_pres_corr.resize(model->wells.size());
    Cf_well.resize(model->wells.size());
    for (int i = 0; i < model->wells.size(); i++)
    {
        pwf_perm_corr[i] = vtkSmartPointer<vtkDoubleArray>::New();
        pwf_pres_corr[i] = vtkSmartPointer<vtkDoubleArray>::New();
        q_perm_corr[i] = vtkSmartPointer<vtkDoubleArray>::New();
        q_pres_corr[i] = vtkSmartPointer<vtkDoubleArray>::New();
        Cf_well[i] = vtkSmartPointer<vtkDoubleArray>::New();
        pwf_perm_corr[i]->SetName(("BHP#" + to_string(i) + " - logperm_corr").c_str());
        pwf_pres_corr[i]->SetName(("BHP#" + to_string(i) + " - pres_corr").c_str());
        q_perm_corr[i]->SetName(("Rate#" + to_string(i) + " - logperm_corr").c_str());
        q_pres_corr[i]->SetName(("Rate#" + to_string(i) + " - pres_corr").c_str());
        Cf_well[i]->SetName(("LogPerm#" + to_string(i) + " - logperm_corr").c_str());
    }

	std::vector<vtkSmartPointer<vtkDoubleArray>> Cp_well;
	Cp_well.resize(model->possible_steps_num);
	for (size_t time_step = 0; time_step < model->possible_steps_num; time_step++)
	{
		Cp_well[time_step] = vtkSmartPointer<vtkDoubleArray>::New();
		Cp_well[time_step]->SetName(("p_corr_well_step#" + to_string(time_step)).c_str());
	}

	auto p_var = vtkSmartPointer<vtkDoubleArray>::New();
	p_var->SetName("pres_variance");
	auto p_std = vtkSmartPointer<vtkDoubleArray>::New();
	p_std->SetName("pres_standart_deviation");
    auto perm = vtkSmartPointer<vtkDoubleArray>::New();
    perm->SetName("Permeability");
    auto perm_var = vtkSmartPointer<vtkDoubleArray>::New();
    perm_var->SetName("Permeability_variance");
    auto perm_stand_dev = vtkSmartPointer<vtkDoubleArray>::New();
    perm_stand_dev->SetName("Permeability_stand_dev");
	auto q_avg_0 = vtkSmartPointer<vtkDoubleArray>::New();
	q_avg_0->SetName("q_avg_0");
	q_avg_0->SetNumberOfComponents(2);
	auto q_avg_2 = vtkSmartPointer<vtkDoubleArray>::New();
	q_avg_2->SetName("q_avg_2");
	q_avg_2->SetNumberOfComponents(2);
	auto qx_std = vtkSmartPointer<vtkDoubleArray>::New();
	qx_std->SetName("qx_standart_deviation");
	auto qy_std = vtkSmartPointer<vtkDoubleArray>::New();
	qy_std->SetName("qy_standart_deviation");

	points->Allocate((num_x + 1) * (num_y + 1));
	cells->Allocate(num_x * num_y);
	double hx = mesh->hx / (double)num_x;
	double hy = mesh->hy / (double)num_y;

	for (int i = 0; i < num_x + 1; i++)
		for (int j = 0; j < num_y + 1; j++)
		{
			const Cell& cell = mesh->cells[(num_y + 2) * i + j];
			points->InsertNextPoint(R_dim * (cell.cent.x + cell.hx / 2), R_dim * (cell.cent.y + cell.hy / 2), 0.0);
		}
	grid->SetPoints(points);

	double q_comps[2];
	size_t x_ind, y_ind;
	double var, Kg, Jx, Jy, dCfp_dx, dCfp_dy, Sigma2;
    double buf1, buf2, buf3, buf4, buf5;
	for (const auto& cell : mesh->cells)
	{
		if (cell.type == elem::QUAD)
		{
			x_ind = cell.id / (num_y + 2) - 1;
			y_ind = cell.id % (num_y + 2) - 1;

			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
			quad->GetPointIds()->SetId(0, y_ind + x_ind * (num_y + 1));
			quad->GetPointIds()->SetId(1, y_ind + x_ind * (num_y + 1) + 1);
			quad->GetPointIds()->SetId(2, y_ind + (x_ind + 1) * (num_y + 1) + 1);
			quad->GetPointIds()->SetId(3, y_ind + (x_ind + 1) * (num_y + 1));
			cells->InsertNextCell(quad);

			p0->InsertNextValue(model->p0_next[cell.id] * model->P_dim / BAR_TO_PA);
			p2->InsertNextValue(model->p2_next[cell.id] * model->P_dim / BAR_TO_PA);
            buf1 = exp(model->getFavg(cell)) * model->props_oil.visc;
            perm->InsertNextValue(M2toMilliDarcy(buf1 * R_dim * R_dim));

			var = model->Cp_next[snap_idx][cell.id * model->cellsNum + cell.id] * model->P_dim / BAR_TO_PA * model->P_dim / BAR_TO_PA;
			p_var->InsertNextValue(var);
			if(var > 0.0)
				p_std->InsertNextValue(sqrt(var));
			else
				p_std->InsertNextValue(0.0);

			const auto& beta_y_minus = mesh->cells[cell.stencil[1]];
			const auto& beta_y_plus = mesh->cells[cell.stencil[2]];
			const auto& beta_x_minus = mesh->cells[cell.stencil[3]];
			const auto& beta_x_plus = mesh->cells[cell.stencil[4]];

			Kg = model->getKg(cell);
			Jx = -(model->p0_next[beta_x_plus.id] - model->p0_next[beta_x_minus.id]) / (beta_x_plus.cent.x - beta_x_minus.cent.x);
			Jy = -(model->p0_next[beta_y_plus.id] - model->p0_next[beta_y_minus.id]) / (beta_y_plus.cent.y - beta_y_minus.cent.y);
			dCfp_dx = (model->Cfp_prev[cell.id * model->cellsNum + beta_x_plus.id] - model->Cfp_prev[cell.id * model->cellsNum + beta_x_minus.id]) / 
						(beta_x_plus.cent.x - beta_x_minus.cent.x);
			dCfp_dy = (model->Cfp_prev[cell.id * model->cellsNum + beta_y_plus.id] - model->Cfp_prev[cell.id * model->cellsNum + beta_y_minus.id]) / 
						(beta_y_plus.cent.y - beta_y_minus.cent.y);
			Sigma2 = model->getSigma2f(cell);

			q_comps[0] = Kg * Jx * cell.hy * cell.hz * model->Q_dim * 86400.0;
			q_comps[1] = Kg * Jy * cell.hx * cell.hz * model->Q_dim * 86400.0;
			q_avg_0->InsertNextTuple(q_comps);
			q_comps[0] = -Kg * ( (model->p2_next[beta_x_plus.id] - model->p2_next[beta_x_minus.id])	/ 
				(beta_x_plus.cent.x - beta_x_minus.cent.x) + dCfp_dx ) * cell.hy * cell.hz * model->Q_dim * 86400.0 - q_comps[0] * Sigma2 / 2.0;
			q_comps[1] = -Kg * ((model->p2_next[beta_y_plus.id] - model->p2_next[beta_y_minus.id]) / 
				(beta_y_plus.cent.y - beta_y_minus.cent.y) + dCfp_dy ) * cell.hx * cell.hz * model->Q_dim * 86400.0 - q_comps[1] * Sigma2 / 2.0;
			q_avg_2->InsertNextTuple(q_comps);

			var = Kg * Kg * (Jx * Jx * Sigma2 - 2.0 * Jx * dCfp_dx + 
			((model->Cp_next[snap_idx][model->cellsNum * beta_x_plus.id + beta_x_plus.id] - model->Cp_next[snap_idx][model->cellsNum * beta_x_minus.id + beta_x_plus.id]) /
				(beta_x_plus.cent.x - beta_x_minus.cent.x) - 
			(model->Cp_next[snap_idx][model->cellsNum * beta_x_plus.id + beta_x_minus.id] - model->Cp_next[snap_idx][model->cellsNum * beta_x_minus.id + beta_x_minus.id]) / 
				(beta_x_plus.cent.x - beta_x_minus.cent.x)) / (beta_x_plus.cent.x - beta_x_minus.cent.x)) * cell.hy * cell.hz * cell.hy * cell.hz * model->Q_dim * 86400.0 * model->Q_dim * 86400.0;
			if (var > 0.0)
				qx_std->InsertNextValue(sqrt(var));
			else
				qx_std->InsertNextValue(0.0);

			var = Kg * Kg * (Jy * Jy * Sigma2 - 2.0 * Jy * dCfp_dy +
			((model->Cp_next[snap_idx][model->cellsNum * beta_y_plus.id + beta_y_plus.id] - model->Cp_next[snap_idx][model->cellsNum * beta_y_minus.id + beta_y_plus.id]) /
				(beta_y_plus.cent.y - beta_y_minus.cent.y) -
			(model->Cp_next[snap_idx][model->cellsNum * beta_y_plus.id + beta_y_minus.id] - model->Cp_next[snap_idx][model->cellsNum * beta_y_minus.id + beta_y_minus.id]) /
				(beta_y_plus.cent.y - beta_y_minus.cent.y)) / (beta_y_plus.cent.y - beta_y_minus.cent.y)) * cell.hx * cell.hz * cell.hx * cell.hz * model->Q_dim * 86400.0 * model->Q_dim * 86400.0;
			if (var > 0.0)
				qy_std->InsertNextValue(sqrt(var));
			else
				qy_std->InsertNextValue(0.0);

            for (size_t time_step = 0; time_step < model->possible_steps_num; time_step++)
            {
                Cp_well[time_step]->InsertNextValue(model->Cp_next[time_step][model->wells.back().cell_id * model->cellsNum + cell.id] /
                                                                        sqrt(model->Cp_next[time_step][model->wells.back().cell_id * model->cellsNum + model->wells.back().cell_id] *
                                                                            model->Cp_next[time_step][cell.id * model->cellsNum + cell.id]));
            }

            buf2 = (exp(model->getSigma2f(cell)) - 1.0) * buf1 * buf1;
            perm_var->InsertNextValue(M2toMilliDarcy(M2toMilliDarcy(buf2 * R_dim * R_dim) * R_dim * R_dim));
            perm_stand_dev->InsertNextValue(M2toMilliDarcy(sqrt(buf2) * R_dim * R_dim));
            for (int i = 0; i < model->wells.size(); i++)
            {
                const auto& well = model->wells[i];
                buf3 = model->Cfp_next[cell.id * model->cellsNum + well.cell_id];
                buf4 = model->getSigma2f(cell);
                buf5 = model->Cp_next[snap_idx][well.cell_id * model->cellsNum + well.cell_id];
                if (fabs(buf3) == 0.0 && (sqrt(buf4) == 0.0 || sqrt(buf5) == 0.0))
                    buf1 = 0.0;
                else
                    buf1 = buf3 / sqrt(buf4 * buf5);
                buf3 = model->Cp_next[snap_idx][well.cell_id * model->cellsNum + cell.id];
                buf4 = model->Cp_next[snap_idx][well.cell_id * model->cellsNum + well.cell_id];
                buf5 = model->Cp_next[snap_idx][cell.id * model->cellsNum + cell.id];
                if (fabs(buf3) == 0.0 && (sqrt(buf4) == 0.0 || sqrt(buf5) == 0.0))
                    buf2 = 0.0;
                else
                    buf2 = buf3 / sqrt(buf4 * buf5);
                if (well.cur_bound)
                {
                    pwf_perm_corr[i]->InsertNextValue(buf1);
                    pwf_pres_corr[i]->InsertNextValue(buf2);
                    q_perm_corr[i]->InsertNextValue(0.0);
                    q_pres_corr[i]->InsertNextValue(0.0);
                }
                else
                {
                    pwf_perm_corr[i]->InsertNextValue(0.0);
                    pwf_pres_corr[i]->InsertNextValue(0.0);
                    q_perm_corr[i]->InsertNextValue(buf1);
                    q_pres_corr[i]->InsertNextValue(buf2);
                }
                Cf_well[i]->InsertNextValue(model->getCf(mesh->cells[well.cell_id], cell) / 
                                    sqrt(model->getSigma2f(mesh->cells[well.cell_id]) * model->getSigma2f(cell)));
            }
		}
	}
	grid->SetCells(VTK_QUAD, cells);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(p0);
	fd->AddArray(p2);
    fd->AddArray(perm);
    for (int i = 0; i < model->wells.size(); i++)
    {
        fd->AddArray(pwf_perm_corr[i]);
        fd->AddArray(pwf_pres_corr[i]);
        fd->AddArray(q_perm_corr[i]);
        fd->AddArray(q_pres_corr[i]);
        fd->AddArray(Cf_well[i]);
    }
	for (const auto& Cp_cur : Cp_well)
		fd->AddArray(Cp_cur);
	fd->AddArray(p_var);
	fd->AddArray(p_std);
	fd->AddArray(q_avg_0);
	fd->AddArray(q_avg_2);
	fd->AddArray(qx_std);
	fd->AddArray(qy_std);
    fd->AddArray(perm_var);
    fd->AddArray(perm_stand_dev);

	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<oil::Oil>;
template class VTKSnapshotter<stoch_oil::StochOil>;
