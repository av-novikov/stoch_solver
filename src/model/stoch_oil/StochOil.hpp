#ifndef STOCH_OIL_HPP_
#define STOCH_OIL_HPP_

#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/stoch_oil/Properties.hpp"
#include "src/Well.hpp"

namespace stoch_oil
{
	typedef var::containers::TapeVar1Phase TapeVariable0;
	typedef var::containers::TapeStochVar1Phase1 TapeVariable1;
	typedef var::containers::TapeStochVar1Phase2 TapeVariable2;
	class StochOil : public AbstractModel<Properties,mesh::RectangularUniformGrid,StochOil,var::StochVariables,
				var::containers::Var1phase,var::containers::StochVar1phase1, var::containers::StochVar1phase2>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class StochOilMethod;
	public:

	protected:
		void makeDimLess();
		void setInitialState();

		static const int var_size0 = VarContainer0::size;
		static const int var_size1 = VarContainer1::size;
		static const int var_size2 = VarContainer2::size;

		TapeVariable0* x0;
		TapeVariable1* x1;
		TapeVariable2* x2;
		adouble* h0;
		adouble* h1;
		adouble* h2;

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
		inline double getS(const Cell& cell) const
		{
			return getPoro(cell) * props_oil.beta + props_sk.beta;
		};

		TapeVariable0 solveInner0(const Cell& cell) const;
		TapeVariable0 solveBorder0(const Cell& cell) const;
		TapeVariable0 solveSource0(const Well& cell) const;

		TapeVariable1 solveInner1(const Cell& cell) const;
		TapeVariable1 solveBorder1(const Cell& cell) const;
		TapeVariable1 solveSource1(const Well& cell) const;
	public:
		StochOil();
		~StochOil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
		double getRate(const Well& well) const;
		double getPwf(const Well& well) const;
	};
};

#endif /* STOCH_OIL_HPP_ */