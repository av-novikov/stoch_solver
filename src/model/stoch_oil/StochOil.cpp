#include "src/model/stoch_oil/StochOil.hpp"

#include <assert.h>

using namespace stoch_oil;

StochOil::StochOil()
{
}
StochOil::~StochOil()
{
}
void StochOil::setProps(const Properties& props)
{
	R_dim = props.R_dim;
	t_dim = props.t_dim;
	Q_dim = R_dim * R_dim * R_dim / t_dim;
	mesh = std::make_shared<Mesh>(*new Mesh(props.num_x, props.num_y, props.hx / R_dim, props.hy / R_dim, props.hz / R_dim));

	cellsNum = mesh.get()->num;
	u_prev0.resize(var_size0 * cellsNum);	u_iter0.resize(var_size0 * cellsNum);	u_next0.resize(var_size0 * cellsNum);
	u_prev1.resize(var_size1 * cellsNum);	u_iter1.resize(var_size1 * cellsNum);	u_next1.resize(var_size1 * cellsNum);
	u_prev2.resize(var_size2 * cellsNum);	u_iter2.resize(var_size2 * cellsNum);	u_next2.resize(var_size2 * cellsNum);
	x0 = new TapeVariable0[cellsNum];		x1 = new TapeVariable1[cellsNum];		x2 = new TapeVariable2[cellsNum];
	h0 = new adouble[var_size0 * cellsNum]; h1 = new adouble[var_size1 * cellsNum]; h2 = new adouble[var_size2 * cellsNum];
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
void StochOil::makeDimLess()
{
	P_dim = props_sk.p_init;

	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	props_sk.p_init /= P_dim;
	props_sk.p_out /= P_dim;
	props_sk.perm /= R_dim * R_dim;
	props_sk.beta /= (1.0 / P_dim);

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
void StochOil::setInitialState()
{
	for (size_t i = 0; i < cellsNum; i++)
	{
		const auto& cell = mesh->cells[i];
		auto data = (*this)[i];
		data.u_prev0.p0 = data.u_iter0.p0 = data.u_next0.p0 = props_sk.p_init;
		
		data.u_prev1.p2 = data.u_iter1.p2 = data.u_next1.p2 = 0.0;
		data.u_prev1.Cfp = data.u_iter1.Cfp = data.u_next1.Cfp = 0.0;

		data.u_prev2.Cp = data.u_iter2.Cp = data.u_next2.Cp = 0.0;
	}

	// WI calculation
	for (auto& well : wells)
	{
		const Cell& cell = mesh->cells[well.cell_id];
		well.r_peaceman = 0.28 * sqrt(cell.hx * cell.hx + cell.hy * cell.hy) / 2.0;
		well.WI = 2.0 * M_PI * props_sk.perm * cell.hz / log(well.r_peaceman / well.rw);
	}
}
void StochOil::setPeriod(const int period)
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
double StochOil::getRate(const Well& well) const
{
	if (well.cur_bound)
		return well.cur_rate;
	else
	{
		const auto& var = (*this)[well.cell_id];
		double p_cell = var.u_next0.p0 + var.u_next0.p0;
		return well.WI * (p_cell - well.cur_pwf) / props_oil.visc;
	}
}
double StochOil::getPwf(const Well& well) const
{
	if (well.cur_bound)
	{
		const auto& var = (*this)[well.cell_id];
		double p_cell = var.u_next0.p0 + var.u_next0.p0;
		return p_cell - well.cur_rate * props_oil.visc / well.WI;
	}
	else
		return well.cur_pwf;
}

TapeVariable0 StochOil::solveInner0(const Cell& cell) const
{
	assert(cell.type == elem::QUAD);

	const auto& next = x0[cell.id];
	const auto prev = (*this)[cell.id].u_prev0;
	
	TapeVariable0 H;
	H.p0 = getS(cell) * (next.p0 - prev.p0) / ht / getKg(cell);

	const auto& beta_y_minus = mesh->cells[cell.stencil[1]];
	const auto& beta_y_plus = mesh->cells[cell.stencil[2]];
	const auto& beta_x_minus = mesh->cells[cell.stencil[3]];
	const auto& beta_x_plus = mesh->cells[cell.stencil[4]];

	const auto& nebr_y_minus = x0[cell.stencil[1]];
	const auto& nebr_y_plus = x0[cell.stencil[2]];
	const auto& nebr_x_minus = x0[cell.stencil[3]];
	const auto& nebr_x_plus = x0[cell.stencil[4]];

	H.p0 -= ((nebr_x_plus.p0 - next.p0) / (beta_x_plus.cent.x - cell.cent.x) - 
			(next.p0 - nebr_x_minus.p0) / (cell.cent.x - beta_x_minus.cent.x)) / cell.hx;

	H.p0 -= (getFavg(beta_x_plus) - getFavg(beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
			(nebr_x_plus.p0 - nebr_x_minus.p0) / (beta_x_plus.cent.x - beta_x_minus.cent.x);

	H.p0 -= ((nebr_y_plus.p0 - next.p0) / (beta_y_plus.cent.y - cell.cent.y) -
		(next.p0 - nebr_y_minus.p0) / (cell.cent.y - beta_y_minus.cent.y)) / cell.hy;

	H.p0 -= (getFavg(beta_y_plus) - getFavg(beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(nebr_y_plus.p0 - nebr_y_minus.p0) / (beta_y_plus.cent.y - beta_y_minus.cent.y);

	return H;
}
TapeVariable0 StochOil::solveBorder0(const Cell& cell) const
{
	assert(cell.type == elem::BORDER);

	TapeVariable0 H;
	const auto& cur = x0[cell.id];
	const auto& nebr = x0[cell.stencil[1]];
	H.p0 = (cur.p0 - (adouble)(props_sk.p_out)) / P_dim;
	return H;
}
TapeVariable0 StochOil::solveSource0(const Well& well) const
{
	const Cell& cell = mesh->cells[well.cell_id];
	TapeVariable0 H;
	H.p0 = well.cur_rate * props_oil.rho_stc / mesh->hz / getKg(cell);
	return H;
}

TapeVariable1 StochOil::solveInner1(const Cell& cell) const
{
	assert(cell.type == elem::QUAD);
	const auto& next = x1[cell.id];
	const auto prev = (*this)[cell.id].u_prev1;

	TapeVariable1 H;
	H.p2 = getS(cell) * (next.p2 - prev.p2) / ht / getKg(cell);

	const auto& beta_y_minus = mesh->cells[cell.stencil[1]];
	const auto& beta_y_plus = mesh->cells[cell.stencil[2]];
	const auto& beta_x_minus = mesh->cells[cell.stencil[3]];
	const auto& beta_x_plus = mesh->cells[cell.stencil[4]];

	const auto& nebr_y_minus = x0[cell.stencil[1]];
	const auto& nebr_y_plus = x0[cell.stencil[2]];
	const auto& nebr_x_minus = x0[cell.stencil[3]];
	const auto& nebr_x_plus = x0[cell.stencil[4]];

	return H;
}
TapeVariable1 StochOil::solveBorder1(const Cell& cell) const
{
	assert(cell.type == elem::BORDER);
	TapeVariable1 H;
	const auto& cur = x1[cell.id];
	H.p2 = cur.p2;
	H.Cfp = cur.Cfp;
	return H;
}
TapeVariable1 StochOil::solveSource1(const Well& cell) const
{
	TapeVariable1 H;
	return H;
}