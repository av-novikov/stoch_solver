#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <math.h>
#include "adolc/adouble.h"

class Interpolate
{
	public:
	Interpolate();
	Interpolate(double *ptx, double *pty, size_t N);
	Interpolate(double *ptx, double *pty, double *dpty, size_t N);
	Interpolate(double *ptx, double *pty, double *dpty,double *d2pty, size_t N);
	~Interpolate();
	
	double Solve(double arg);
	adouble Solve(adouble arg);
	double DSolve(double arg);
	double D2Solve(double arg);
	
	double *x;
	double *y;
	double *dy;
	double *d2y;
	int *Flag;
	size_t N;
	size_t Nf;
	double xmin, xmax, temp;
	bool IsSetDiff;
	bool IsSetDiff2;
};

#endif /* INTERPOLATE_H_ */