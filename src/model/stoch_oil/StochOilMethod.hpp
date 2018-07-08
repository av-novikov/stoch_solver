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
		void writeData();

		void solveStep();
		void solveStep_p0();
		void solveStep_Cfp();
		void solveStep_p2();
		void solveStep_Cp();

		std::ofstream plot_P, plot_Q, pvd;
		ParSolver solver0, solver1;
		int step_idx;
		double averVal, averValPrev, dAverVal;

		static const int var_size = 1;
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

		void computeJac_p0();
		void computeJac_Cfp(const int cell_id);
		void computeJac_p2();
		void computeJac_Cp(const int cell_id);
		void fillIndices0();
		void fillIndices1();
		void fill_p0();
		void fill_Cfp(const int cell_id);
		void fill_p2();
		void fill_Cp(const int cell_id);
		void copySolution_p0(const paralution::LocalVector<double>& sol);
		void copySolution_Cfp(const int cell_id, const paralution::LocalVector<double>& sol);
		void copySolution_p2(const paralution::LocalVector<double>& sol);
		void copySolution_Cp(const int cell_id, const paralution::LocalVector<double>& sol);

		void copyTimeLayer();
		void copyIterLayer_p0();
		void copyIterLayer_p2();
		void copyIterLayer_Cp(const int cell_id);

		double convergance_p0(int& ind, int& varInd);
		double convergance_p2(int& ind, int& varInd);
		double convergance_Cp(int& ind, int& varInd, const int cell_id);

		double averValue_p0() const;
		double averValue_Cfp(const int cell_id) const;
		double averValue_p2() const;
		double averValue_Cp(const int cell_id) const;
	public:
		StochOilMethod(Model* _model);
		~StochOilMethod();

		void start();
	};
};

#endif /* STOCH_OILMETHOD_HPP_ */