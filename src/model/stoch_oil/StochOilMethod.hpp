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

		static const int var_size = 1;
		double** jac_p0;
		double* y_p0;
		int* ind_i_p0;
		int* ind_j_p0;
		double* a_p0;
		int* ind_rhs_p0;
		double* rhs_p0;
		int* cols_p0;
		// Number of non-zero elements in sparse matrix
		int elemNum_p0;

		double** jac_Cfp;
		double* y_Cfp;
		int* ind_i_Cfp;
		int* ind_j_Cfp;
		double* a_Cfp;
		int* ind_rhs_Cfp;
		double* rhs_Cfp;
		int* cols_Cfp;
		// Number of non-zero elements in sparse matrix
		int elemNum_Cfp;

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