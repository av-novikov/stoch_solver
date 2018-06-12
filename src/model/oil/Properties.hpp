#ifndef OIL_PROPERTIES_HPP_
#define OIL_PROPERTIES_HPP_

#include <vector>
#include <utility>
#include "src/Well.hpp"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace oil
{
	struct Skeleton_Props
	{
		double p_init;
		double perm;
	};
	struct Oil_Props
	{
		double visc;
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

#endif /* OIL_PROPERTIES_HPP_ */