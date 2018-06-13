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

	double** jac;
	double* y;

	Model* model;
	Mesh* mesh;

	size_t curTimePeriod;
	const double Tt;
	double cur_t, cur_t_log;

	void copyIterLayer();
	void copyTimeLayer();

	double convergance(int& ind, int& varInd);
	void averValue(std::array<double, var_size>& aver);

	virtual void writeData() = 0;
	virtual void control() = 0;
	virtual void doNextStep();
	virtual void solveStep() = 0;

	std::vector<int> stencil_idx;
	inline void getMatrixStencil(const Cell& cell)
	{
		if (cell.type == elem::BORDER)
		{
			stencil_idx.resize(2);
			stencil_idx[0] = cell.id;
			const size_t ind_x = int(cell.id / (mesh->num_y + 2));
			const size_t ind_y = cell.id % (mesh->num_y + 2);
			if(ind_y == 0)
				stencil_idx[1] = cell.id + 1;
			else if(ind_y == mesh->num_y + 1)
				stencil_idx[1] = cell.id - 1;
			if (ind_x == 0)
				stencil_idx[1] = cell.id + mesh->num_y + 2;
			else if (ind_x == mesh->num_x + 1)
				stencil_idx[1] = cell.id - mesh->num_y - 2;
		}
		else
		{
			stencil_idx.resize(Mesh::stencil);
			stencil_idx[0] = cell.id;
			stencil_idx[1] = cell.id - 1;
			stencil_idx[2] = cell.id + 1;
			stencil_idx[3] = cell.id - mesh->num_y - 2;
			stencil_idx[4] = cell.id + mesh->num_y + 2;
		}
	};

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

	virtual void fill() = 0;
	void fillIndices()
	{
		int counter = 0;

		for (int i = 0; i < model->varNum; i++)
		{
			const auto& cell = mesh->cells[i];
			getMatrixStencil(cell);
			for (size_t i = 0; i < var_size; i++)
				for (const int idx : stencil_idx)
					for (size_t j = 0; j < var_size; j++)
					{
						ind_i[counter] = var_size * cell.id + i;			ind_j[counter++] = var_size * idx + j;
					}
			stencil_idx.clear();
		}

		elemNum = counter;

		for (int i = 0; i < var_size * model->cellsNum; i++)
			ind_rhs[i] = i;
	};
	virtual void start();
};

#endif /* ABSTRACTMETHOD_HPP_ */