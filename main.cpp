#include <iostream>

#include "src/Scene.hpp"
#include "src/model/oil/OilMethod.hpp"
#include "src/model/stoch_oil/StochOilMethod.hpp"

using namespace std;

namespace issues
{
	template <class modelType, class methodType>
	struct Issue
	{
		typename typedef modelType Model;
		typename typedef methodType Method;
	};

	struct Oil : public Issue<oil::Oil, oil::OilMethod> {};
	struct StochOil : public Issue<stoch_oil::StochOil, stoch_oil::StochOilMethod> {};
}

int main()
{
	stoch_oil::Properties props;
	
	props.possible_steps_num = 20;
	props.t_dim = 3600.0;
	props.ht = props.ht_min = 1000.0;
	props.ht_max = 1000000.0;

	props.hx = props.hy = props.R_dim = 2100.0;		props.hz = 10.0;
	props.num_x = props.num_y = 21;
	size_t num = (props.num_x + 2) * (props.num_y + 2);

	props.props_sk.p_init = props.props_sk.p_out = 100.0 * BAR_TO_PA;
	props.props_sk.perm = 50.0;
	props.props_sk.m = 0.1;
	props.props_sk.beta = 4.E-10;
	props.props_sk.l_f = 100.0;
	props.props_sk.sigma_f = 0.5;

	props.props_oil.visc = 1.0;
	props.props_oil.rho_stc = 887.261;
	props.props_oil.beta = 1.0 * 1.e-9;
	props.props_oil.p_ref = props.props_sk.p_init;

	props.wells.push_back(Well(0, (props.num_y + 2) * int(3 * (props.num_x + 2) / 4) + 3 * (props.num_y + 2) / 4));
	auto& well1 = props.wells.back();
	well1.periodsNum = 1;
	well1.period.resize(well1.periodsNum);
	well1.period[0] = 31.0 * 86400.0;
	well1.rate.resize(well1.periodsNum);
	well1.rate[0] = 50.0;
	well1.leftBoundIsRate.resize(well1.periodsNum);
	well1.leftBoundIsRate[0] = true;
	well1.rw = 0.1;

	props.wells.push_back(Well(1, (props.num_y + 2) * int((props.num_x + 2) / 4) + (props.num_y + 2) / 4));
	auto& well2 = props.wells.back();
	well2.periodsNum = 1;
	well2.period.resize(well2.periodsNum);
	well2.period[0] = 31.0 * 86400.0;
	well2.rate.resize(well2.periodsNum);
	well2.rate[0] = 50.0;
	well2.leftBoundIsRate.resize(well2.periodsNum);
	well2.leftBoundIsRate[0] = true;
	well2.rw = 0.1;

	Scene<issues::StochOil> scene;
	scene.load(props);
	scene.start();

	return 0;
}