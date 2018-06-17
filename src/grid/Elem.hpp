#ifndef ELEM_HPP_
#define ELEM_HPP_

#include "src/grid/Point.hpp"

namespace elem
{
	typedef point::Point Point;

	enum Type {QUAD, BORDER};
	template <int MaxStencilNum>
	class Element
	{
	protected:
		
	public:
		int id;
		const Type type;
		static const int getCurStencNum()
		{
			if (type == QUAD)
				return MaxStencilNum;
			else if (type == BORDER)
				return 2;
			else
				return -1;
		};
		std::array<int, MaxStencilNum> stencil;
		std::array<double, MaxStencilNum - 1> trans;

		Element(const Type _type) : type(_type) {};
		Element(const int _id, const Type _type) : id(_id), type(_type) {};
	};
	class Quad : public Element<5>
	{
	public:
		Point cent;

		double hx, hy, hz;
		double V;

		Quad(const int _id, const Type _type) : Element(_id, _type) {};
		Quad(const int _id, const Type _type, const Point _cent, const Point _sizes) : 
			Element(_id, _type), cent(_cent), hx(_sizes.x), hy(_sizes.y), hz(_sizes.z), V(_sizes.x * _sizes.y * _sizes.z) {};
	};
};

#endif /* ELEM_HPP_ */
