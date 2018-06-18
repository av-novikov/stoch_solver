#ifndef OIL_HPP_
#define OIL_HPP_

#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/oil/Properties.hpp"
#include "src/Well.hpp"

namespace oil
{
	typedef var::containers::TapeVar1Phase TapeVariable;
	class Oil : public AbstractModel<Properties, mesh::RectangularUniformGrid,Oil,var::BasicVariables,var::containers::Var1phase>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class OilMethod;
	protected:
		void makeDimLess();
		void setInitialState();

		TapeVariable* x;
		adouble* h;

		Skeleton_Props props_sk;
		Oil_Props props_oil;
		std::vector<Well> wells;

		inline double getPoro(const Cell& cell) const
		{
			return props_sk.m;
		};
		inline double getPerm(const Cell& cell) const
		{
			return props_sk.perm;
		};
		inline double getFavg(const Cell& cell) const
		{
			return log(getPerm(cell) / props_oil.visc);
		};
		inline double getKg(const Cell& cell) const
		{
			return exp(getFavg(cell));
		};

		adouble solveInner(const Cell& cell) const;
		adouble solveBorder(const Cell& cell) const;
	public:
		Oil();
		~Oil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
		double getRate(const Well& well) const;
		double getPwf(const Well& well) const;
	};
};

#endif /* OIL_HPP_ */