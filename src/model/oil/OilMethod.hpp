#ifndef OILMETHOD_HPP_
#define OILMETHOD_HPP_

#include "src/model/AbstractMethod.hpp"
#include "src/model/oil/Oil.hpp"
#include "src/utils/ParalutionInterface.h"

namespace oil
{
	class OilMethod : public AbstractMethod<Oil>
	{
	protected:
		void control();
		void solveStep();
		void writeData();

		std::ofstream plot_P, plot_Q;
		ParSolver solver;

		void computeJac();
		void fill();
		void copySolution(const paralution::LocalVector<double>& sol);
	public:
		OilMethod(Model* _model);
		~OilMethod();

		void start();
	};
};

#endif /* OILMETHOD_HPP_ */