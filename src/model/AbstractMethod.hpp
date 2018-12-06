#ifndef ABSTRACTMETHOD_HPP_
#define ABSTRACTMETHOD_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <array>
#include <vector>

template <class modelType>
class AbstractMethod {
public:
	typedef modelType Model;
	typedef typename Model::Mesh Mesh;
	typedef typename Model::Cell Cell;
	static const int var_size = Model::var_size;
protected:
	double t_dim;
	const size_t size;

	Model* model;
	Mesh* mesh;

	size_t curTimePeriod;
	const double Tt;
	double cur_t, cur_t_log;

	virtual void copyIterLayer();
	virtual void copyTimeLayer();

	double convergance(int& ind, int& varInd);
	void averValue(std::array<double, var_size>& aver);

	virtual void writeData() = 0;
	virtual void control() = 0;
	virtual void doNextStep();
	virtual void solveStep() = 0;

	inline void getMatrixStencil(Cell& cell)
	{
		const size_t ind_x = int(cell.id / (mesh->num_y + 2));
		const size_t ind_y = cell.id % (mesh->num_y + 2);
		cell.stencil[0] = cell.id;

		const double k1 = model->getPerm(cell);
		double k2;

		if (cell.type == elem::BORDER)
		{
			if (ind_y == 0)
			{
				cell.stencil[1] = cell.id + 1;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
			}
			else if (ind_y == mesh->num_y + 1)
			{
				cell.stencil[1] = cell.id - 1;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
			}
			if (ind_x == 0)
			{
				cell.stencil[1] = cell.id + mesh->num_y + 2;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
			}
			else if (ind_x == mesh->num_x + 1)
			{
				cell.stencil[1] = cell.id - mesh->num_y - 2;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
			}
		}
		else
		{
			cell.stencil[1] = cell.id - 1;
			const Cell& beta1 = mesh->cells[cell.stencil[1]];
			k2 = model->getPerm(beta1);
			cell.trans[0] = k1 * k2 * (cell.hy + beta1.hy) / (k1 * beta1.hy + k2 * cell.hy);
			
            cell.stencil[2] = cell.id + 1;
			const Cell& beta2 = mesh->cells[cell.stencil[2]];
			k2 = model->getPerm(beta2);
			cell.trans[1] = k1 * k2 * (cell.hy + beta2.hy) / (k1 * beta2.hy + k2 * cell.hy);
			
            cell.stencil[3] = cell.id - mesh->num_y - 2;
			const Cell& beta3 = mesh->cells[cell.stencil[3]];
			k2 = model->getPerm(beta3);
			cell.trans[2] = k1 * k2 * (cell.hx + beta3.hx) / (k1 * beta3.hx + k2 * cell.hx);
			
            cell.stencil[4] = cell.id + mesh->num_y + 2;
			const Cell& beta4 = mesh->cells[cell.stencil[4]];
			k2 = model->getPerm(beta4);
			cell.trans[3] = k1 * k2 * (cell.hx + beta4.hx) / (k1 * beta4.hx + k2 * cell.hx);
		}
	};

	double** jac;
	double* y;
	int* ind_i;
	int* ind_j;
	double* a;
	int* ind_rhs;
	double* rhs;
	int* cols;
	// Number of non-zero elements in sparse matrix
	int elemNum;

	int options[4];
	int repeat;

public:
	AbstractMethod(modelType* _model);
	virtual ~AbstractMethod();

	virtual void fill() {};
	void fillIndices()
	{
		int counter = 0;

		for (int i = 0; i < model->varNum; i++)
		{
			auto& cell = mesh->cells[i];
			getMatrixStencil(cell);
			
			for (size_t i = 0; i < var_size; i++)
				for (const int idx : cell.stencil)
					for (size_t j = 0; j < var_size; j++)
					{
						ind_i[counter] = var_size * cell.id + i;			ind_j[counter++] = var_size * idx + j;
					}
		}

		elemNum = counter;

		for (int i = 0; i < var_size * model->cellsNum; i++)
			ind_rhs[i] = i;
	};
	virtual void start();
};

#endif /* ABSTRACTMETHOD_HPP_ */