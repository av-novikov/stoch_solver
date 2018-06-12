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
	virtual void fillIndices() = 0;
	virtual void start();
};

#endif /* ABSTRACTMETHOD_HPP_ */