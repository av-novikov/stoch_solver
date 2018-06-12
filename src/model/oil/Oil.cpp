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
	t_dim = props.t_dim;
	Q_dim = R_dim * R_dim * R_dim / t_dim;
	mesh = std::make_shared<Mesh>(*new Mesh(props.num_x, props.num_y, props.hx / R_dim, props.hy / R_dim, props.hz / R_dim));

	cellsNum = mesh.get()->num;
	varNum = VarContainer::size * cellsNum;
	u_prev.resize(varNum);
	u_iter.resize(varNum);
	u_next.resize(varNum);
	x = new TapeVariable[cellsNum];
	h = new adouble[var_size * cellsNum];
	Volume = mesh.get()->V;

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_sk = props.props_sk;
	props_sk.perm = MilliDarcyToM2(props_sk.perm);

	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	wells = props.wells;

	makeDimLess();
}
void Oil::makeDimLess()
{
	P_dim = props_sk.p_init;

	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	props_sk.p_init /= P_dim;
	props_sk.perm /= R_dim * R_dim;

	for (auto& well : wells)
	{
		well.period /= t_dim;
		well.pwf /= P_dim;
		well.rate /= Q_dim;
	}

	props_oil.visc /= P_dim * t_dim;
}
void Oil::setInitialState()
{
	for (size_t i = 0; i < cellsNum; i++)
	{
		const auto& cell = mesh->cells[i];
		auto data = (*this)[i];
		data.u_prev.p0 = data.u_iter.p0 = data.u_next.p0 = props_sk.p_init;
	}
}
void Oil::setPeriod(const int period)
{
	for (auto& well : wells)
	{
		well.cur_bound = well.leftBoundIsRate[period];
		if (well.cur_bound)
			well.cur_rate = well.rate[period];
		else
			well.cur_pwf = well.pwf[period];
	}
}