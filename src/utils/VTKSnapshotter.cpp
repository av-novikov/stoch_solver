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
	
	auto Cfp_well = vtkSmartPointer<vtkDoubleArray>::New();
	Cfp_well->SetName("Cfp_well");
	
	/*std::vector<vtkSmartPointer<vtkDoubleArray>> Cfp;
	Cfp.resize(mesh->num_x * mesh->num_y);
	int counter = 0;
	for (const auto& cell : mesh->cells)
		if (cell.type == elem::QUAD)
		{
			Cfp[counter] = vtkSmartPointer<vtkDoubleArray>::New();
			Cfp[counter++]->SetName(("Cfp#" + to_string(cell.id)).c_str());
		}*/

	std::vector<vtkSmartPointer<vtkDoubleArray>> Cp_well;
	Cp_well.resize(model->possible_steps_num);
	for (size_t time_step = 0; time_step < model->possible_steps_num; time_step++)
	{
		Cp_well[time_step] = vtkSmartPointer<vtkDoubleArray>::New();
		Cp_well[time_step]->SetName(("Cp_well #" + to_string(time_step)).c_str());
	}

	auto variance = vtkSmartPointer<vtkDoubleArray>::New();
	variance->SetName("variance");
	auto stand_dev = vtkSmartPointer<vtkDoubleArray>::New();
	stand_dev->SetName("standart_deviation");

	auto q_avg_0 = vtkSmartPointer<vtkDoubleArray>::New();
	q_avg_0->SetName("q_avg_0");
	q_avg_0->SetNumberOfComponents(2);
	auto q_avg_2 = vtkSmartPointer<vtkDoubleArray>::New();
	q_avg_2->SetName("q_avg_2");
	q_avg_2->SetNumberOfComponents(2);

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
	double var;
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
			Cfp_well->InsertNextValue(model->Cfp_next[model->wells.back().cell_id * model->cellsNum + cell.id] * model->P_dim / BAR_TO_PA);

			var = model->Cp_next[snap_idx][cell.id * model->cellsNum + cell.id] * model->P_dim / BAR_TO_PA * model->P_dim / BAR_TO_PA;
			variance->InsertNextValue(var);
			if(var > 0.0)
				stand_dev->InsertNextValue(sqrt(var));
			else
				stand_dev->InsertNextValue(0.0);

			const auto& beta_y_minus = mesh->cells[cell.stencil[1]];
			const auto& beta_y_plus = mesh->cells[cell.stencil[2]];
			const auto& beta_x_minus = mesh->cells[cell.stencil[3]];
			const auto& beta_x_plus = mesh->cells[cell.stencil[4]];
			q_comps[0] = -model->getKg(cell) * (model->p0_next[beta_x_plus.id] - model->p0_next[beta_x_minus.id]) / 
												(beta_x_plus.cent.x - beta_x_minus.cent.x);
			q_comps[1] = -model->getKg(cell) * (model->p0_next[beta_y_plus.id] - model->p0_next[beta_y_minus.id]) /
				(beta_y_plus.cent.y - beta_y_minus.cent.y);
			q_avg_0->InsertNextTuple(q_comps);
			q_comps[0] = -model->getKg(cell) * ( (model->p2_next[beta_x_plus.id] - model->p2_next[beta_x_minus.id] + 
				model->Cfp_prev[cell.id * model->cellsNum + beta_x_plus.id] - model->Cfp_prev[cell.id * model->cellsNum + beta_x_minus.id]) 
				/ (beta_x_plus.cent.x - beta_x_minus.cent.x) ) - q_comps[0] * model->getSigma2f(cell) / 2.0;
			q_comps[1] = -model->getKg(cell) * ((model->p2_next[beta_y_plus.id] - model->p2_next[beta_y_minus.id] +
				model->Cfp_prev[cell.id * model->cellsNum + beta_y_plus.id] - model->Cfp_prev[cell.id * model->cellsNum + beta_y_minus.id])
				/ (beta_y_plus.cent.y - beta_y_minus.cent.y)) - q_comps[1] * model->getSigma2f(cell) / 2.0;
			q_avg_2->InsertNextTuple(q_comps);

			/*counter = 0;
			for (const auto& cur_cell : mesh->cells)
				if (cur_cell.type == elem::QUAD)
				{
					auto& Cfp_cur = Cfp[counter++];
					Cfp_cur->InsertNextValue(model->Cfp_prev[cur_cell.id * model->cellsNum + cell.id]);
				}*/

			for (size_t time_step = 0; time_step < model->possible_steps_num; time_step++)
				Cp_well[time_step]->InsertNextValue(model->Cp_next[time_step][model->wells.back().cell_id * model->cellsNum + cell.id]);
		}
	}
	grid->SetCells(VTK_QUAD, cells);

	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(p0);
	fd->AddArray(p2);
	fd->AddArray(Cfp_well);
	//for(const auto& Cfp_cur : Cfp)
	//	fd->AddArray(Cfp_cur);
	for (const auto& Cp_cur : Cp_well)
		fd->AddArray(Cp_cur);
	fd->AddArray(variance);
	fd->AddArray(stand_dev);
	fd->AddArray(q_avg_0);
	fd->AddArray(q_avg_2);

	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(snap_idx).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<oil::Oil>;
template class VTKSnapshotter<stoch_oil::StochOil>;
