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
	mesh = std::make_shared<grid::Mesh>(new Mesh(props.num_x, props.num_y, props.hx / R_dim, props.hy / R_dim, props.hz / R_dim));
}
void Oil::loadMesh()
{
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