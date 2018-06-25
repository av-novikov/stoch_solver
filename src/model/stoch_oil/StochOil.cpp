#include "src/model/stoch_oil/StochOil.hpp"

#include <valarray>

#include <assert.h>

using namespace stoch_oil;
using std::valarray;

StochOil::StochOil()
{
}
StochOil::~StochOil()
{
	delete[] x, h;
}
void StochOil::setProps(const Properties& props)
{
	R_dim = props.R_dim;
	t_dim = props.t_dim;
	Q_dim = R_dim * R_dim * R_dim / t_dim;
	mesh = std::make_shared<Mesh>(*new Mesh(props.num_x, props.num_y, props.hx / R_dim, props.hy / R_dim, props.hz / R_dim));

	cellsNum = mesh.get()->num;
	p0_prev.resize(var_size * cellsNum);	
	p0_iter.resize(var_size * cellsNum);	
	p0_next.resize(var_size * cellsNum);

	Cfp_prev.resize(var_size * cellsNum * cellsNum);	
	Cfp_iter.resize(var_size * cellsNum * cellsNum);
	Cfp_next.resize(var_size * cellsNum * cellsNum);

	p2_prev.resize(var_size * cellsNum);	
	p2_iter.resize(var_size * cellsNum);	
	p2_next.resize(var_size * cellsNum);

	Cp_prev.resize(var_size * cellsNum * cellsNum);	
	Cp_iter.resize(var_size * cellsNum * cellsNum);	
	Cp_next.resize(var_size * cellsNum * cellsNum);

	x = new adouble[cellsNum];		
	h = new adouble[cellsNum]; 
	Volume = mesh.get()->V;

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_sk = props.props_sk;
	props_sk.perm = MilliDarcyToM2(props_sk.perm);

	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	wells = props.wells;
	for (auto& well : wells)
		for (auto& rate : well.rate)
			rate /= 86400;

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
	props_sk.l_f /= R_dim;

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
		p0_prev[i] = p0_iter[i] = p0_next[i] = props_sk.p_init;
		p2_prev[i] = p2_iter[i] = p2_next[i] = 0.0;

		const auto sl = std::slice(i * cellsNum, cellsNum, var_size);
		for (size_t j = i * cellsNum; j < (i + 1) * cellsNum; j++)
		{
			Cfp_prev[j] = Cfp_iter[j] = Cfp_next[j] = 0.0;
			Cp_prev[j] = Cp_iter[j] = Cp_next[j] = 0.0;
		}
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
		well.cur_period = period;
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
		double p_cell = p0_next[well.cell_id];
		return well.WI * (p_cell - well.cur_pwf) / props_oil.visc;
	}
}
double StochOil::getPwf(const Well& well) const
{
	if (well.cur_bound)
	{
		double p_cell = p0_next[well.cell_id];
		return p_cell - well.cur_rate * props_oil.visc / well.WI;
	}
	else
		return well.cur_pwf;
}

adouble StochOil::solveInner_p0(const Cell& cell) const
{
	assert(cell.type == elem::QUAD);

	const auto& next = x[cell.id];
	const auto prev = p0_prev[cell.id];
	
	adouble H;
	H = getS(cell) * (next - prev) / getKg(cell);

	const auto& beta_y_minus = mesh->cells[cell.stencil[1]];
	const auto& beta_y_plus = mesh->cells[cell.stencil[2]];
	const auto& beta_x_minus = mesh->cells[cell.stencil[3]];
	const auto& beta_x_plus = mesh->cells[cell.stencil[4]];

	const auto& nebr_y_minus = x[cell.stencil[1]];
	const auto& nebr_y_plus = x[cell.stencil[2]];
	const auto& nebr_x_minus = x[cell.stencil[3]];
	const auto& nebr_x_plus = x[cell.stencil[4]];

	H -= ht * ((nebr_x_plus - next) / (beta_x_plus.cent.x - cell.cent.x) -
			(next - nebr_x_minus) / (cell.cent.x - beta_x_minus.cent.x)) / cell.hx;

	H -= ht * (getFavg(beta_x_plus) - getFavg(beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
			(nebr_x_plus - nebr_x_minus) / (beta_x_plus.cent.x - beta_x_minus.cent.x);

	H -= ht * ((nebr_y_plus - next) / (beta_y_plus.cent.y - cell.cent.y) -
		(next - nebr_y_minus) / (cell.cent.y - beta_y_minus.cent.y)) / cell.hy;

	H -= ht * (getFavg(beta_y_plus) - getFavg(beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(nebr_y_plus - nebr_y_minus) / (beta_y_plus.cent.y - beta_y_minus.cent.y);

	return H;
}
adouble StochOil::solveBorder_p0(const Cell& cell) const
{
	assert(cell.type == elem::BORDER);

	const auto& cur = x[cell.id];
	const auto& nebr = x[cell.stencil[1]];
	return (cur - (adouble)(props_sk.p_out)) / P_dim;
}
adouble StochOil::solveSource_p0(const Well& well) const
{
	const Cell& cell = mesh->cells[well.cell_id];
	return well.cur_rate * ht / cell.V / getKg(cell);
}

adouble StochOil::solveInner_Cfp(const Cell& cell, const Cell& cur_cell) const
{
	assert(cell.type == elem::QUAD && cur_cell.type == elem::QUAD);
	const auto& next = x[cell.id];
	const auto prev = Cfp_prev[cur_cell.id * cellsNum + cell.id];

	adouble H;
	H = getS(cell) * (next - prev) / getKg(cell);

	const int& y_minus = cell.stencil[1];
	const int& y_plus = cell.stencil[2];
	const int& x_minus = cell.stencil[3];
	const int& x_plus = cell.stencil[4];

	const auto& beta_y_minus = mesh->cells[y_minus];
	const auto& beta_y_plus = mesh->cells[y_plus];
	const auto& beta_x_minus = mesh->cells[x_minus];
	const auto& beta_x_plus = mesh->cells[x_plus];

	const auto& nebr_y_minus = x[y_minus];
	const auto& nebr_y_plus = x[y_plus];
	const auto& nebr_x_minus = x[x_minus];
	const auto& nebr_x_plus = x[x_plus];

	H -= ht * ((nebr_x_plus - next) / (beta_x_plus.cent.x - cell.cent.x) -
		(next - nebr_x_minus) / (cell.cent.x - beta_x_minus.cent.x)) / cell.hx;

	H -= ht * (getFavg(beta_x_plus) - getFavg(beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
		(nebr_x_plus - nebr_x_minus) / (beta_x_plus.cent.x - beta_x_minus.cent.x);

	H -= ht * ((nebr_y_plus - next) / (beta_y_plus.cent.y - cell.cent.y) -
		(next - nebr_y_minus) / (cell.cent.y - beta_y_minus.cent.y)) / cell.hy;

	H -= ht * (getFavg(beta_y_plus) - getFavg(beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(nebr_y_plus - nebr_y_minus) / (beta_y_plus.cent.y - beta_y_minus.cent.y);


	double H1 = -ht * ((p0_next[x_plus] - p0_next[x_minus]) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
	(getCf(cur_cell, beta_x_plus) - getCf(cur_cell, beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) +
					(p0_next[y_plus] - p0_next[y_minus]) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
	(getCf(cur_cell, beta_y_plus) - getCf(cur_cell, beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y));
	
	double H2 = -getS(cell) / getKg(cell) * (p0_next[cell.id] - p0_prev[cell.id]) * getCf(cur_cell, cell);

	return H + H1 + H2;
}
adouble StochOil::solveBorder_Cfp(const Cell& cell, const Cell& cur_cell) const
{
	assert(cell.type == elem::BORDER || cur_cell.type == elem::BORDER);
	return x[cell.id] / P_dim;
}
adouble StochOil::solveSource_Cfp(const Well& well, const Cell& cur_cell) const
{
	const Cell& cell = mesh->cells[well.cell_id];
	return -well.cur_rate * ht / cell.V / getKg(cell) * getCf(cur_cell, cell);
}

adouble StochOil::solveInner_p2(const Cell& cell) const
{
	assert(cell.type == elem::QUAD);

	const auto& next = x[cell.id];
	const auto prev = p2_prev[cell.id];

	adouble H;
	H = getS(cell) * (next - prev) / getKg(cell);

	const int& y_minus = cell.stencil[1];
	const int& y_plus = cell.stencil[2];
	const int& x_minus = cell.stencil[3];
	const int& x_plus = cell.stencil[4];

	const auto& beta_y_minus = mesh->cells[y_minus];
	const auto& beta_y_plus = mesh->cells[y_plus];
	const auto& beta_x_minus = mesh->cells[x_minus];
	const auto& beta_x_plus = mesh->cells[x_plus];

	const auto& nebr_y_minus = x[y_minus];
	const auto& nebr_y_plus = x[y_plus];
	const auto& nebr_x_minus = x[x_minus];
	const auto& nebr_x_plus = x[x_plus];

	H -= ht * ((nebr_x_plus - next) / (beta_x_plus.cent.x - cell.cent.x) -
		(next - nebr_x_minus) / (cell.cent.x - beta_x_minus.cent.x)) / cell.hx;

	H -= ht * (getFavg(beta_x_plus) - getFavg(beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
		(nebr_x_plus - nebr_x_minus) / (beta_x_plus.cent.x - beta_x_minus.cent.x);

	H -= ht * ((nebr_y_plus - next) / (beta_y_plus.cent.y - cell.cent.y) -
		(next - nebr_y_minus) / (cell.cent.y - beta_y_minus.cent.y)) / cell.hy;

	H -= ht * (getFavg(beta_y_plus) - getFavg(beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(nebr_y_plus - nebr_y_minus) / (beta_y_plus.cent.y - beta_y_minus.cent.y);

	double H1 = -ht * ((Cfp_next[(cell.id + 1) * cellsNum + x_plus] - Cfp_next[(cell.id + 1) * cellsNum + x_minus]) -
						(Cfp_next[(cell.id - 1) * cellsNum + x_plus] - Cfp_next[(cell.id - 1) * cellsNum + x_minus])) /
						(beta_x_plus.cent.x - beta_x_minus.cent.x) / (beta_x_plus.cent.x - beta_x_minus.cent.x) - 
				ht * ((Cfp_next[(cell.id + 1) * cellsNum + y_plus] - Cfp_next[(cell.id + 1) * cellsNum + y_minus]) -
						(Cfp_next[(cell.id - 1) * cellsNum + y_plus] - Cfp_next[(cell.id - 1) * cellsNum + y_minus])) /
						(beta_y_plus.cent.y - beta_y_minus.cent.y) / (beta_y_plus.cent.y - beta_y_minus.cent.y);

	const size_t idx = cell.id * cellsNum + cell.id;
	double H2 = getS(cell) / getKg(cell) * ((p0_next[cell.id] - p0_prev[cell.id]) * getSigmaf(cell) / 2.0 -
							(Cfp_next[idx] - Cfp_prev[idx]));
	return H + H1 + H2;
}
adouble StochOil::solveBorder_p2(const Cell& cell) const
{
	assert(cell.type == elem::BORDER);
	return x[cell.id] / P_dim;
}
adouble StochOil::solveSource_p2(const Well& well) const
{
	const Cell& cell = mesh->cells[well.cell_id];
	return well.cur_rate * ht / cell.V / getKg(cell) * getSigmaf(cell) / 2.0;
}

adouble StochOil::solveInner_Cp(const Cell& cell, const Cell& cur_cell) const
{
	assert(cell.type == elem::QUAD && cur_cell.type == elem::QUAD);
	const auto& next = x[cell.id];
	const auto prev = Cp_prev[cur_cell.id * cellsNum + cell.id];

	adouble H;
	H = getS(cell) * (next - prev) / getKg(cell);

	const int& y_minus = cell.stencil[1];
	const int& y_plus = cell.stencil[2];
	const int& x_minus = cell.stencil[3];
	const int& x_plus = cell.stencil[4];

	const auto& beta_y_minus = mesh->cells[y_minus];
	const auto& beta_y_plus = mesh->cells[y_plus];
	const auto& beta_x_minus = mesh->cells[x_minus];
	const auto& beta_x_plus = mesh->cells[x_plus];

	const auto& nebr_y_minus = x[y_minus];
	const auto& nebr_y_plus = x[y_plus];
	const auto& nebr_x_minus = x[x_minus];
	const auto& nebr_x_plus = x[x_plus];

	H -= ht * ((nebr_x_plus - next) / (beta_x_plus.cent.x - cell.cent.x) -
		(next - nebr_x_minus) / (cell.cent.x - beta_x_minus.cent.x)) / cell.hx;

	H -= ht * (getFavg(beta_x_plus) - getFavg(beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
		(nebr_x_plus - nebr_x_minus) / (beta_x_plus.cent.x - beta_x_minus.cent.x);

	H -= ht * ((nebr_y_plus - next) / (beta_y_plus.cent.y - cell.cent.y) -
		(next - nebr_y_minus) / (cell.cent.y - beta_y_minus.cent.y)) / cell.hy;

	H -= ht * (getFavg(beta_y_plus) - getFavg(beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(nebr_y_plus - nebr_y_minus) / (beta_y_plus.cent.y - beta_y_minus.cent.y);


	double H1 = -ht * ((p0_next[x_plus] - p0_next[x_minus]) / (beta_x_plus.cent.x - beta_x_minus.cent.x) *
		(Cfp_next[]getCf(cur_cell, beta_x_plus) - getCf(cur_cell, beta_x_minus)) / (beta_x_plus.cent.x - beta_x_minus.cent.x) +
		(p0_next[y_plus] - p0_next[y_minus]) / (beta_y_plus.cent.y - beta_y_minus.cent.y) *
		(getCf(cur_cell, beta_y_plus) - getCf(cur_cell, beta_y_minus)) / (beta_y_plus.cent.y - beta_y_minus.cent.y));

	double H2 = -getS(cell) / getKg(cell) * (p0_next[cell.id] - p0_prev[cell.id]) * getCf(cur_cell, cell);

	return H + H1 + H2;
}
adouble StochOil::solveBorder_Cp(const Cell& cell, const Cell& cur_cell) const
{
	assert(cell.type == elem::BORDER || cur_cell.type == elem::BORDER);
	return x[cell.id] / P_dim;
}
adouble StochOil::solveSource_Cp(const Well& well, const Cell& cur_cell) const
{
	const Cell& cell = mesh->cells[well.cell_id];
	return -well.cur_rate * ht / cell.V / getKg(cell) * Cfp_next[cur_cell.id * cellsNum + cell.id];
}

