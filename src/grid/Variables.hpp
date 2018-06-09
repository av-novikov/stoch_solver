#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <valarray>
#include <array>
#include "adolc/adouble.h"

namespace var
{
	namespace containers
	{
		struct Var1phase
		{
			static const int size = 1;
			double& p0;

			Var1phase(double* data) : p0(data[0]) {};
			Var1phase(const double* data) : p0(const_cast<double&>(data[0])) {};
		};
		struct TapeVar1Phase
		{
			static const int size = 1;
			adouble p0;
		};
	};

	template <typename TVariable>
	struct BaseVarWrapper
	{
		TVariable u_prev, u_iter, u_next;
	};
	template <typename TVariable>
	struct BasicVariables
	{
		static const int size = TVariable::size;
		typedef BaseVarWrapper<TVariable> Wrap;
		std::valarray<double> u_prev, u_iter, u_next;

		Wrap operator[](const size_t idx)
		{
			return{ TVariable(&u_prev[idx * size]), TVariable(&u_iter[idx * size]), TVariable(&u_next[idx * size]) };
		};
		Wrap operator[](const size_t idx) const
		{
			return{ TVariable(&u_prev[idx * size]), TVariable(&u_iter[idx * size]), TVariable(&u_next[idx * size]) };
		};
	};
}

#endif /* VARIABLES_HPP_ */