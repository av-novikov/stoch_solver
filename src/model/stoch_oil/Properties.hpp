#ifndef STOCH_OIL_PROPERTIES_HPP_
#define STOCH_OIL_PROPERTIES_HPP_

#include <vector>
#include <utility>
#include "src/Well.hpp"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace stoch_oil
{
	struct Skeleton_Props
	{
		double m;
		double p_init;
		double p_out;
		double perm;
		double beta;
	};
	struct Oil_Props
	{
		double visc;
		double rho_stc;
		double beta;

		double p_ref;
		inline adouble getB(adouble p) const
		{
			return exp(-(adouble)beta * (p - p_ref));
		};
		inline adouble getDensity(adouble p) const
		{
			return rho_stc / getB(p);
		};
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};

	};
	struct Properties
	{
		double ht, ht_min, ht_max;

		Skeleton_Props props_sk;
		Oil_Props props_oil;
		std::vector<Well> wells;

		double R_dim, t_dim;

		int num_x, num_y;
		double hx, hy, hz;
	};
};

#endif /* STOCH_OIL_PROPERTIES_HPP_ */