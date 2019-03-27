#ifndef MESH_HPP_
#define MESH_HPP_

#include <src/grid/Elem.hpp>

namespace mesh
{
	typedef elem::Point Point;

	class CellRectangularUniformGrid
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
	public:
		typedef elem::Quad Cell;
		static const int stencil = 5;
	public:
		std::vector<Cell> cells;
		
		const int num_x, num_y, num;
		const double hx, hy, hz;
		const double V;
	public:
		CellRectangularUniformGrid(const int _num_x, const int _num_y, const double _hx, const double _hy, const double _hz) :
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
		~CellRectangularUniformGrid() {};
	};
    class NodeRectangularUniformGrid
    {
        template<typename> friend class snapshotter::VTKSnapshotter;
        template<typename> friend class AbstractMethod;
    public:
        typedef elem::DualQuad Node;
        static const int stencil = 5;

        std::vector<Node> nodes;
        const int num_x, num_y, num;
        const double hx, hy, hz;
        const double V;

        NodeRectangularUniformGrid(const int _num_x, const int _num_y, const double _hx, const double _hy, const double _hz) :
            num_x(_num_x), num_y(_num_y), num((_num_x + 1) * (_num_y + 1)), hx(_hx), hy(_hy), hz(_hz), V(_hx * _hy * _hz)
        {
            Point cur(0.0, 0.0, hz / 2.0);
            double hx1 = hx / (double)num_x;
            double hy1 = hy / (double)num_y;
            const Point size(hx1, hy1, hz);
            Point cur_size;
            int counter = 0;
            elem::Type ntype;
            for (int i = 0; i < num_x + 1; i++)
            {
                for (int j = 0; j < num_y + 1; j++)
                {
                    if ((i == 0 && j == 0) || (i == 0 && j == num_y) || (i == num_x && j == 0) || (i == num_x && j == num_y))
                    {
                        cur_size = { hx1 / 2.0, hy1 / 2.0, hz };
                        ntype = elem::CORNER;
                    }
                    else if (i == 0 || i == num_x)
                    {
                        cur_size = { hx1 / 2.0, hy1, hz };
                        ntype = elem::BORDER;
                    } 
                    else if (j == 0 || j == num_y)
                    {
                        cur_size = { hx1, hy1 / 2.0, hz };
                        ntype = elem::BORDER;
                    }
                    else
                    {
                        cur_size = { hx1, hy1, hz };
                        ntype = elem::QUAD;
                    }

                    nodes.push_back(Node(counter++, ntype, cur, cur_size));
                    cur.y += hy1;
                }
                cur.x += hx1;
            }
        };
        ~NodeRectangularUniformGrid() {};
    };
};

#endif /* MESH_HPP_ */