#ifndef WELL_HPP_
#define WELL_HPP_

#define _USE_MATH_DEFINES
#include <math.h>

#include <valarray>

class Well
{
public:
	size_t cell_id;
	// Number of periods
	size_t periodsNum;
	// End times of periods [sec]
	std::valarray<double> period;
	// Oil rates [m3/day]
	std::valarray<double> rate;
	// Vector of BHPs [bar]
	std::valarray<double> pwf;
	// If left boundary condition would be 2nd type
	std::valarray<bool> leftBoundIsRate;

	size_t cur_period;
	double cur_rate, cur_pwf;
	bool cur_bound;

	double rw, r_peaceman, WI;

	Well(const size_t _cell_id) : cell_id(_cell_id) {};
	~Well() {};
};

#endif /* WELL_HPP_ */
