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

		const double k1 = model->getGeomPerm(cell);
		double k2;

		if (cell.type == elem::BORDER)
		{
			if (ind_y == 0)
			{
				cell.stencil[1] = cell.id + 1;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getGeomPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
			}
			else if (ind_y == mesh->num_y + 1)
			{
				cell.stencil[1] = cell.id - 1;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getGeomPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
			}
			if (ind_x == 0)
			{
				cell.stencil[1] = cell.id + mesh->num_y + 2;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getGeomPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
			}
			else if (ind_x == mesh->num_x + 1)
			{
				cell.stencil[1] = cell.id - mesh->num_y - 2;
				const Cell& beta = mesh->cells[cell.stencil[1]];
				k2 = model->getGeomPerm(beta);
				cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
			}
		}
		else
		{
			cell.stencil[1] = cell.id - 1;
			const Cell& beta1 = mesh->cells[cell.stencil[1]];
			k2 = model->getGeomPerm(beta1);
			cell.trans[0] = k1 * k2 * (cell.hy + beta1.hy) / (k1 * beta1.hy + k2 * cell.hy);
			
            cell.stencil[2] = cell.id + 1;
			const Cell& beta2 = mesh->cells[cell.stencil[2]];
			k2 = model->getGeomPerm(beta2);
			cell.trans[1] = k1 * k2 * (cell.hy + beta2.hy) / (k1 * beta2.hy + k2 * cell.hy);
			
            cell.stencil[3] = cell.id - mesh->num_y - 2;
			const Cell& beta3 = mesh->cells[cell.stencil[3]];
			k2 = model->getGeomPerm(beta3);
			cell.trans[2] = k1 * k2 * (cell.hx + beta3.hx) / (k1 * beta3.hx + k2 * cell.hx);
			
            cell.stencil[4] = cell.id + mesh->num_y + 2;
			const Cell& beta4 = mesh->cells[cell.stencil[4]];
			k2 = model->getGeomPerm(beta4);
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

template <class modelType>
class AbstractDualGridMethod {
public:
    typedef modelType Model;
    typedef typename Model::CellMesh CellMesh;
    typedef typename Model::NodeMesh NodeMesh;
    typedef typename Model::Cell Cell;
    typedef typename Model::Node Node;
    static const int var_size = Model::var_size;
protected:
    double t_dim;
    const size_t cells_size, nodes_size;

    Model* model;
    CellMesh* cell_mesh;
    NodeMesh* node_mesh;

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

    inline void getCellMatrixStencil(Cell& cell)
    {
        const size_t ind_x = int(cell.id / (cell_mesh->num_y + 2));
        const size_t ind_y = cell.id % (cell_mesh->num_y + 2);
        cell.stencil[0] = cell.id;

        const double k1 = model->getGeomPerm(cell);
        double k2;

        if (cell.type == elem::BORDER)
        {
            if (ind_y == 0)
            {
                cell.stencil[1] = cell.id + 1;
                const Cell& beta = cell_mesh->cells[cell.stencil[1]];
                k2 = model->getGeomPerm(beta);
                cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
            }
            else if (ind_y == cell_mesh->num_y + 1)
            {
                cell.stencil[1] = cell.id - 1;
                const Cell& beta = cell_mesh->cells[cell.stencil[1]];
                k2 = model->getGeomPerm(beta);
                cell.trans[0] = k1 * k2 * (cell.hy + beta.hy) / (k1 * beta.hy + k2 * cell.hy);
            }
            if (ind_x == 0)
            {
                cell.stencil[1] = cell.id + cell_mesh->num_y + 2;
                const Cell& beta = cell_mesh->cells[cell.stencil[1]];
                k2 = model->getGeomPerm(beta);
                cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
            }
            else if (ind_x == cell_mesh->num_x + 1)
            {
                cell.stencil[1] = cell.id - cell_mesh->num_y - 2;
                const Cell& beta = cell_mesh->cells[cell.stencil[1]];
                k2 = model->getGeomPerm(beta);
                cell.trans[0] = k1 * k2 * (cell.hx + beta.hx) / (k1 * beta.hx + k2 * cell.hx);
            }
        }
        else
        {
            cell.stencil[1] = cell.id - 1;
            const Cell& beta1 = cell_mesh->cells[cell.stencil[1]];
            k2 = model->getGeomPerm(beta1);
            cell.trans[0] = k1 * k2 * (cell.hy + beta1.hy) / (k1 * beta1.hy + k2 * cell.hy);

            cell.stencil[2] = cell.id + 1;
            const Cell& beta2 = cell_mesh->cells[cell.stencil[2]];
            k2 = model->getGeomPerm(beta2);
            cell.trans[1] = k1 * k2 * (cell.hy + beta2.hy) / (k1 * beta2.hy + k2 * cell.hy);

            cell.stencil[3] = cell.id - cell_mesh->num_y - 2;
            const Cell& beta3 = cell_mesh->cells[cell.stencil[3]];
            k2 = model->getGeomPerm(beta3);
            cell.trans[2] = k1 * k2 * (cell.hx + beta3.hx) / (k1 * beta3.hx + k2 * cell.hx);

            cell.stencil[4] = cell.id + cell_mesh->num_y + 2;
            const Cell& beta4 = cell_mesh->cells[cell.stencil[4]];
            k2 = model->getGeomPerm(beta4);
            cell.trans[3] = k1 * k2 * (cell.hx + beta4.hx) / (k1 * beta4.hx + k2 * cell.hx);
        }
    };
    inline void getNodeMatrixStencil(Node& node)
    {
        const auto& mesh = node_mesh;
        const size_t ind_x = int(node.id / (mesh->num_y + 1));
        const size_t ind_y = node.id % (mesh->num_y + 1);
        node.stencil[0] = node.id;
        double k1, k2;

        if (node.type == elem::CORNER)
        { 
            if (ind_y == 0 && ind_x == 0)
            {
                node.stencil[1] = node.id + 1;
                node.stencil[2] = node.id + mesh->num_y + 1;
                const Cell& beta = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y + 1];
                k1 = model->getGeomPerm(beta);
                node.trans[0] = node.trans[1] = k1;
            }
            else if (ind_y == mesh->num_y && ind_x == 0)
            {
                node.stencil[1] = node.id - 1;
                node.stencil[2] = node.id + mesh->num_y + 1;
                const Cell& beta = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y];
                node.trans[0] = node.trans[1] = model->getGeomPerm(beta);
            }
            else if (ind_y == 0 && ind_x == mesh->num_x)
            {
                node.stencil[1] = node.id + 1;
                node.stencil[2] = node.id - mesh->num_y - 1;
                const Cell& beta = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y + 1];
                node.trans[0] = node.trans[1] = model->getGeomPerm(beta);
            }
            else if (ind_y == mesh->num_y && ind_x == mesh->num_x)
            {
                node.stencil[1] = node.id - 1;
                node.stencil[2] = node.id - mesh->num_y - 1;
                const Cell& beta = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y];
                node.trans[0] = node.trans[1] = model->getGeomPerm(beta);
            }
        }
        else if (node.type == elem::BORDER)
        {

            if (ind_y == 0)
            {
                node.stencil[1] = node.id + 1;
                const Cell& beta1 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y + 1];
                const Cell& beta2 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y + 1];
                k1 = model->getGeomPerm(beta1);     k2 = model->getGeomPerm(beta2);
                node.trans[0] = (k1 * beta1.hx + k2 * beta2.hx) / (beta1.hx + beta2.hx);
            }
            else if (ind_y == mesh->num_y)
            {
                node.stencil[1] = node.id - 1;
                const Cell& beta1 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y];
                const Cell& beta2 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y];
                k1 = model->getGeomPerm(beta1);     k2 = model->getGeomPerm(beta2);
                node.trans[0] = (k1 * beta1.hx + k2 * beta2.hx) / (beta1.hx + beta2.hx);
            }
            if (ind_x == 0)
            {
                node.stencil[1] = node.id + mesh->num_y + 1;
                const Cell& beta1 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y];
                const Cell& beta2 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y + 1];
                k1 = model->getGeomPerm(beta1);     k2 = model->getGeomPerm(beta2);
                node.trans[0] = (k1 * beta1.hy + k2 * beta2.hy) / (beta1.hy + beta2.hy);
            }
            else if (ind_x == mesh->num_x)
            {
                node.stencil[1] = node.id - mesh->num_y - 1;
                const Cell& beta1 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y];
                const Cell& beta2 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y + 1];
                k1 = model->getGeomPerm(beta1);     k2 = model->getGeomPerm(beta2);
                node.trans[0] = (k1 * beta1.hy + k2 * beta2.hy) / (beta1.hy + beta2.hy);
            }
        }
        else if (node.type == elem::QUAD)
        {
            node.stencil[1] = node.id - 1;
            const Cell& beta11 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y];
            const Cell& beta12 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y];
            k1 = model->getGeomPerm(beta11);     k2 = model->getGeomPerm(beta12);
            node.trans[0] = (k1 * beta11.hx + k2 * beta12.hx) / (beta11.hx + beta12.hx);

            node.stencil[2] = node.id + 1;
            const Cell& beta21 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y + 1];
            const Cell& beta22 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y + 1];
            k1 = model->getGeomPerm(beta21);     k2 = model->getGeomPerm(beta22);
            node.trans[1] = (k1 * beta21.hx + k2 * beta22.hx) / (beta21.hx + beta22.hx);

            node.stencil[3] = node.id - mesh->num_y - 1;
            const Cell& beta31 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y];
            const Cell& beta32 = cell_mesh->cells[ind_x * (cell_mesh->num_y + 2) + ind_y + 1];
            k1 = model->getGeomPerm(beta31);     k2 = model->getGeomPerm(beta32);
            node.trans[2] = (k1 * beta31.hy + k2 * beta32.hy) / (beta31.hy + beta32.hy);

            node.stencil[4] = node.id + mesh->num_y + 1;
            const Cell& beta41 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y];
            const Cell& beta42 = cell_mesh->cells[(ind_x + 1) * (cell_mesh->num_y + 2) + ind_y + 1];
            k1 = model->getGeomPerm(beta41);     k2 = model->getGeomPerm(beta42);
            node.trans[3] = (k1 * beta41.hy + k2 * beta42.hy) / (beta41.hy + beta42.hy);
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
    AbstractDualGridMethod(modelType* _model);
    virtual ~AbstractDualGridMethod();

    virtual void fill() {};
    virtual void start();
};


#endif /* ABSTRACTMETHOD_HPP_ */