#ifndef OIL_PROPERTIES_HPP_
#define OIL_PROPERTIES_HPP_

#include <vector>
#include <utility>

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace oil
{
	struct Skeleton_Props
	{
	};
	struct Oil_Props
	{
	};
	struct Properties
	{
		double R_dim;

		int num_x, num_y;
		double hx, hy, hz;
	};
};

#endif /* OIL_PROPERTIES_HPP_ */