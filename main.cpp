#include <iostream>

#include "src/Scene.hpp"
#include "src/model/oil/OilMethod.hpp"

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
}

int main()
{
	oil::Properties props;
	
	props.t_dim = 3600.0;
	props.hx = props.hy = props.R_dim = 100.0;		props.hz = 1.0;
	props.num_x = props.num_y = 10;
	size_t num = (props.num_x + 2) * (props.num_y + 2);

	props.props_sk.p_init = 100.0 * BAR_TO_PA;
	props.props_sk.perm = 100.0;

	props.props_oil.visc = 1.0;

	props.wells.push_back(Well(num / 2));
	auto& well = props.wells.back();
	well.periodsNum = 1;
	well.period.resize(well.periodsNum);
	well.period[0] = 100.0 * 3600.0;
	well.pwf.resize(well.periodsNum);
	well.pwf[0] = 90.0 * BAR_TO_PA;
	well.leftBoundIsRate.resize(well.periodsNum);
	well.leftBoundIsRate[0] = false;

	Scene<issues::Oil> scene;
	scene.load(props);
	scene.start();

	return 0;
}