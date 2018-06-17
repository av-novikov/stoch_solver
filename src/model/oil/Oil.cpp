#include "src/model/oil/Oil.hpp"

#include <assert.h>

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
	for (int i = 0; i < wells.back().periodsNum; i++)
	{
		auto& well = wells[i];
		if (well.leftBoundIsRate[i])
			well.rate[i] /= 86400.0;
	}

	makeDimLess();
}
void Oil::makeDimLess()
{
	P_dim = props_sk.p_init;

	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	props_sk.p_init /= P_dim;
	props_sk.p_out /= P_dim;
	props_sk.perm /= R_dim * R_dim;

	for (auto& well : wells)
	{
		well.rw /= R_dim;
		well.period /= t_dim;
		well.pwf /= P_dim;
		well.rate /= Q_dim;
	}

	props_oil.visc /= P_dim * t_dim;
	props_oil.p_ref /= P_dim;
	props_oil.beta /= (1.0 / P_dim);
	props_oil.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}
void Oil::setInitialState()
{
	for (size_t i = 0; i < cellsNum; i++)
	{
		const auto& cell = mesh->cells[i];
		auto data = (*this)[i];
		data.u_prev.p0 = data.u_iter.p0 = data.u_next.p0 = props_sk.p_init;
	}

	// WI calculation
	for (auto& well : wells)
	{
		const Cell& cell = mesh->cells[well.cell_id];
		well.r_peaceman = 0.28 * sqrt(cell.hx * cell.hx + cell.hy * cell.hy) / 2.0;
		well.WI = 2.0 * M_PI * props_sk.perm * cell.hz / log(well.r_peaceman / well.rw);
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
double Oil::getRate(const Well& well) const
{
	if (well.cur_bound)
		return well.cur_rate;
	else
	{
		double p_cell = (*this)[well.cell_id].u_next.p0;
		return well.WI * (p_cell - well.cur_pwf) / props_oil.visc;
	}
}
double Oil::getPwf(const Well& well) const
{
	if (well.cur_bound)
	{
		double p_cell = (*this)[well.cell_id].u_next.p0;
		return p_cell - well.cur_rate * props_oil.visc / well.WI;
	}
	else
		return well.cur_pwf;
}

adouble Oil::solveInner(const Cell& cell) const
{
	assert(cell.type == elem::QUAD);

	const auto& next = x[cell.id];
	const auto prev = (*this)[cell.id].u_prev;
	adouble H = getPoro(cell) * props_oil.getDensity(next.p0) - getPoro(cell) * props_oil.getDensity(prev.p0);
	adouble interp_middle;

	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = mesh->cells[cell.stencil[i + 1]];
		const auto& nebr = x[cell.stencil[i + 1]];
		if (abs(cell.id - beta.id) == 1)
			interp_middle = linearInterp1d(props_oil.getDensity(next.p0) / props_oil.getViscosity(next.p0), cell.hy,
				props_oil.getDensity(nebr.p0) / props_oil.getViscosity(nebr.p0), beta.hy);
		else
			interp_middle = linearInterp1d(props_oil.getDensity(next.p0) / props_oil.getViscosity(next.p0), cell.hx,
				props_oil.getDensity(nebr.p0) / props_oil.getViscosity(nebr.p0), beta.hx);
		double inter = interp_middle.value();
		H += ht / cell.V * cell.trans[i] * interp_middle * (next.p0 - nebr.p0);
	}

	return H;
}
adouble Oil::solveBorder(const Cell& cell) const
{
	assert(cell.type == elem::BORDER);
	const auto& cur = x[cell.id];
	const auto& nebr = x[cell.stencil[1]];
	return (cur.p0 - (adouble)(props_sk.p_out)) / P_dim;
}