#ifndef RECTANGULAR_UNIFORM_GRID_HPP_
#define RECTANGULAR_UNIFORM_GRID_HPP_

#include <src/grid/Elem.hpp>

namespace mesh
{
	typedef elem::Point Point;

	class RectangularUniformGrid
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
	public:
		typedef elem::Quad Cell;
	public:
		std::vector<Cell> cells;
		
		const int num_x, num_y, num;
		const double hx, hy, hz;
		const double V;
	public:
		RectangularUniformGrid(const int _num_x, const int _num_y, const double _hx, const double _hy, const double _hz) :
			num_x(_num_x), num_y(_num_y), num((_num_x + 2) * (_num_y + 2)), hx(_hx), hy(_hy), hz(_hz), V(_hx * _hy * _hz) 
		{
			Point cur(0.0, 0.0, hz / 2.0);
			double hx1 = hx / (double)num_x;
			double hy1 = hy / (double)num_y;
			const Point size(hx1, hy1, hz);
			int counter = 0;

			cur.y = 0.0;
			cells.push_back(Cell(counter++, elem::BORDER, cur, {0.0, 0.0, hz}));
			cur.y += hy1 / 2;
			for (int j = 1; j < num_y + 1; j++)
			{
				cells.push_back(Cell(counter++, elem::BORDER, cur, {0.0, hy1, hz}));
				cur.y += hy1;
			}
			cur.y -= hy1 / 2;
			cells.push_back(Cell(counter++, elem::BORDER, cur, {0.0, 0.0, hz}));

			cur.x += hx1 / 2;
			for (int i = 1; i < num_x + 1; i++)
			{
				cur.y = 0.0;
				cells.push_back(Cell(counter++, elem::BORDER, cur, {hx1, 0.0, hz}));
				cur.y += hy1 / 2;
				for (int j = 1; j < num_y + 1; j++)
				{
					cells.push_back(Cell(counter++, elem::QUAD, cur, size));
					cur.y += hy1;
				}
				cur.y -= hy1 / 2;
				cells.push_back(Cell(counter++, elem::BORDER, cur, {hx1, 0.0, hz}));

				cur.x += hx1;
			}

			cur.x -= hx1 / 2;
			cur.y = 0.0;
			cells.push_back(Cell(counter++, elem::BORDER, cur, { 0.0, 0.0, hz }));
			cur.y += hy1 / 2;
			for (int j = 1; j < num_y + 1; j++)
			{
				cells.push_back(Cell(counter++, elem::BORDER, cur, { 0.0, hy1, hz }));
				cur.y += hy1;
			}
			cur.y -= hy1 / 2;
			cells.push_back(Cell(counter++, elem::BORDER, cur, { 0.0, 0.0, hz }));

		};
		~RectangularUniformGrid() {};
	};
}

#endif /* RECTANGULAR_UNIFORM_GRID_HPP_ */