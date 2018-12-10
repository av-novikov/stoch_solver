#ifndef WELL_HPP_
#define WELL_HPP_

#define _USE_MATH_DEFINES
#include <math.h>

#include <valarray>

class Well
{
public:
	int id;
	int cell_id;
	// Number of periods
	int periodsNum;
	// End times of periods [sec]
	std::valarray<double> period;
	// Oil rates [m3/day]
	std::valarray<double> rate;
	// Vector of BHPs [bar]
	std::valarray<double> pwf;
	// If left boundary condition would be 2nd type
	std::valarray<bool> leftBoundIsRate;

	int cur_period;
	double cur_rate, cur_pwf;
	bool cur_bound;

	double rw, r_peaceman, WI, perm;

	Well(const int _id, const int _cell_id) : id(_id), cell_id(_cell_id) {};
	~Well() {};
};

#endif /* WELL_HPP_ */
