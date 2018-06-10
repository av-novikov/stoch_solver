#include "src/model/oil/Oil.hpp"

using namespace oil;

Oil::Oil()
{
}
Oil::~Oil()
{
}
void Oil::setProps(const Properties& props)
{
	R_dim = props.R_dim;
	mesh = std::make_shared<Mesh>(*new Mesh(props.num_x, props.num_y, props.hx / R_dim, props.hy / R_dim, props.hz / R_dim));

	/*mshreader::MshReader reader;
	mesh = std::make_shared<grid::Mesh>(*reader.read(nebrFileName, R_dim));
	mesh->process_geometry();

	cellsNum = mesh.get()->getCalcCellsSize();
	varNum = VarContainer::size * cellsNum;
	u_prev.resize(varNum);
	u_iter.resize(varNum);
	u_next.resize(varNum);
	Volume = mesh->Volume;*/
}
void Oil::makeDimLess()
{
}
void Oil::setInitialState()
{
}
void Oil::setPeriod(const int period)
{
}