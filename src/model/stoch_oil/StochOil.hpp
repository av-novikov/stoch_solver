#ifndef STOCH_OIL_HPP_
#define STOCH_OIL_HPP_

#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/stoch_oil/Properties.hpp"
#include "src/Well.hpp"

namespace stoch_oil
{
	/*typedef var::containers::TapeVar1Phase TapeVariable0;
	typedef var::containers::TapeStochVar1Phase1 TapeVariable1;
	typedef var::containers::TapeStochVar1Phase2 TapeVariable2;
	typedef var::containers::TapeStochVar1Phase3 TapeVariable3;*/
	class StochOil : public AbstractModel<Properties,mesh::RectangularUniformGrid,StochOil,var::StochVariables,
				var::containers::Var1phase>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class StochOilMethod;
	public:

	protected:
		void makeDimLess();
		void setInitialState();

		adouble* x_p0;
		adouble* x_Cfp;
		adouble* x_p2;
		adouble* x_Cp;

		adouble* h_p0;
		adouble* h_Cfp;
		adouble* h_p2;
		adouble* h_Cp;

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
			return getPoro(cell) * (props_oil.beta + props_sk.beta);
		};
		inline double getCf(const Cell& cell, const Cell& beta) const
		{
			auto corrFoo1 = [this](const auto& p1, const auto& p2) -> double
			{
				return props_sk.sigma_f * props_sk.sigma_f * exp(-point::distance(p1, p2) / props_sk.l_f);
			};
			auto corrFoo2 = [this](const auto& p1, const auto& p2) -> double
			{
				const double dist = point::distance(p1, p2) / props_sk.l_f;
				return props_sk.sigma_f * props_sk.sigma_f * exp(-dist * dist);
			};
			return corrFoo1(cell.cent, beta.cent);
		};

		adouble solveInner0(const Cell& cell) const;
		adouble solveBorder0(const Cell& cell) const;
		adouble solveSource0(const Well& well) const;

		adouble solveInner1(const Cell& cur_cell, const Cell& cell) const;
		adouble solveBorder1(const Cell& cur_cell, const Cell& cell) const;
		adouble solveSource1(const Well& well, const Cell& cur_cell) const;
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