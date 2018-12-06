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
	
	props.possible_steps_num = 14;
	props.start_time_simple_approx = 4;
	props.t_dim = 3600.0;
	props.ht = props.ht_min = 5000.0;
	props.ht_max = 1000000.0;

	props.hx = props.hy = props.R_dim = 2100.0;		props.hz = 10.0;
	props.num_x = props.num_y = 11;
	size_t num = (props.num_x + 2) * (props.num_y + 2);

	props.props_sk.p_init = props.props_sk.p_out = 159.48 * BAR_TO_PA;
	props.props_sk.perm = 100.0;
	props.props_sk.m = 0.1;
	props.props_sk.beta = 4.E-10;
	props.props_sk.l_f = 500.0;
	props.props_sk.sigma_f = 0.5;

	props.props_oil.visc = 1.0;
	props.props_oil.rho_stc = 887.261;
	props.props_oil.beta = 1.0 * 1.e-9;
	props.props_oil.p_ref = props.props_sk.p_init;

	props.wells.push_back(Well(0, (props.num_y + 2) * 6 + 6));
	auto& well1 = props.wells.back();
	well1.periodsNum = 1;
	well1.period.resize(well1.periodsNum);
	well1.period[0] = 365.0 * 86400.0;
	well1.rate.resize(well1.periodsNum);
	well1.rate[0] = 150.0;
	well1.pwf.resize(well1.periodsNum);
	//well1.pwf[0] = 82.0 * BAR_TO_PA;
	well1.leftBoundIsRate.resize(well1.periodsNum);
	well1.leftBoundIsRate[0] = true;
	well1.rw = 0.1;

    props.conditions.push_back({ well1.cell_id, props.props_sk.perm });
    //props.conditions.push_back({ 3 + 3 * (props.num_x + 2), 1.2 * props.props_sk.perm });

	/*props.wells.push_back(Well(1, (props.num_y + 2) * (props.num_x + 1 - 4) + (props.num_y + 1 - 4)));
	auto& well2 = props.wells.back();
	well2.periodsNum = 1;
	well2.period.resize(well2.periodsNum);
	well2.period[0] = 365.0 * 86400.0;
	well2.rate.resize(well2.periodsNum);
	well2.rate[0] = 150.0;
	well2.leftBoundIsRate.resize(well2.periodsNum);
	well2.leftBoundIsRate[0] = true;
	well2.rw = 0.1;*/

	Scene<issues::StochOil> scene;
	scene.load(props);
	scene.start();

	return 0;
}