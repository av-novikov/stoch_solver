#ifndef RECTANGULAR_UNIFORM_GRID_HPP_
#define RECTANGULAR_UNIFORM_GRID_HPP_

#include <src/grid/Elem.hpp>

namespace mesh
{
	typedef elem::Point Point;

	class RectangularUniformGrid
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractSolver;
	public:
		typedef elem::Quad Cell;
		std::vector<Cell> cells;
		
		const int num_x, num_y, num;
		const double hx, hy, hz;
		const double V;

		RectangularUniformGrid(const int _num_x, const int _num_y, const double _hx, const double _hy, const double _hz) :
			num_x(_num_x), num_y(_num_y), num((_num_x + 2) * (_num_y + 2)), hx(_hx), hy(_hy), hz(_hz), V(_hx * _hy * _hz) 
		{
			Point cur(0.0, 0.0, hz / 2.0);
			const Point size(hx / (double)num_x, hy / (double)num_y, hz);
			elem::Type type;
			int counter = 0;
			for (int i = 0; i < num_x + 2; i++)
			{
				cur.y = 0.0;
				if (i == 1 || i == num_x + 1)
					cur.x += size.x / 2;
				else
					cur.x += size.x;

				for (int j = 0; j < num_y + 2; j++)
				{
					if (j == 1 || j == num_y + 1)
					{
						cur.y += size.y / 2;
						type = elem::BORDER;
					}
					else
					{
						cur.y += size.y;
						type = elem::QUAD;
					}
					if (i == 0 || i == num_x + 1)
						type = elem::BORDER;

					cells.push_back(Cell(counter++, type, cur, size));
				}
			}
		};
		~RectangularUniformGrid() {};
	};
}

#endif /* RECTANGULAR_UNIFORM_GRID_HPP_ */