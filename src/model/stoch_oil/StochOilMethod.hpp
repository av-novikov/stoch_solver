#ifndef STOCH_OILMETHOD_HPP_
#define STOCH_OILMETHOD_HPP_

#include "src/model/AbstractMethod.hpp"
#include "src/model/stoch_oil/StochOil.hpp"
#include "src/utils/ParalutionInterface.h"

namespace stoch_oil
{
	class StochOilMethod : public AbstractMethod<StochOil>
	{
	protected:
		void control();
		void solveStep();
		void writeData();

		std::ofstream plot_P, plot_Q, pvd;
		ParSolver solver0, solver1;
		int step_idx;
		std::array<double, var_size> averVal, averValPrev, dAverVal;

		static const int var_size0 = 1;
		double** jac0;
		double* y0;
		int* ind_i0;
		int* ind_j0;
		double* a0;
		int* ind_rhs0;
		double* rhs0;
		int* cols0;
		// Number of non-zero elements in sparse matrix
		int elemNum0;
		static const int var_size1 = 2;
		double** jac1;
		double* y1;
		int* ind_i1;
		int* ind_j1;
		double* a1;
		int* ind_rhs1;
		double* rhs1;
		int* cols1;
		// Number of non-zero elements in sparse matrix
		int elemNum1;

		void computeJac0();
		void computeJac1();
		//void computeJac1();
		void fillIndices0();
		void fillIndices1();
		void fill0();
		void fill1();
		void copySolution0(const paralution::LocalVector<double>& sol);
		void copySolution1(const paralution::LocalVector<double>& sol);
	public:
		StochOilMethod(Model* _model);
		~StochOilMethod();

		void start();
	};
};

#endif /* STOCH_OILMETHOD_HPP_ */