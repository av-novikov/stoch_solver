#include <iostream>

#include "src/Scene.hpp"
#include "src/model/oil/OilMethod.hpp"
#include "src/model/stoch_oil/StochOilMethod.hpp"
#include "src/model/dual_stoch_oil/DualStochOilMethod.hpp"

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
    struct DualStochOil : public Issue<dual_stoch_oil::DualStochOil, dual_stoch_oil::DualStochOilMethod> {};
}


void loadWells(const double x1, const double x2, const double y1, const double y2, 
                const int num_x, const int num_y, const std::string fileName, 
                std::vector<Well>& wells, std::vector<dual_stoch_oil::Measurement>& conds, double visc)
{
    const double hx = (x2 - x1) / (double)num_x;
    const double hy = (y2 - y1) / (double)num_y;

    std::ifstream file(fileName.c_str(), std::ifstream::in);
    std::string buf;
    std::string::size_type sz;
    int id;
    double x, y, perm;
    int id_x, id_y;
    file >> buf;
    while (!file.eof())
    {
        id = (int)(std::stod(buf, &sz));
        file >> buf;    x = std::stod(buf, &sz);
        file >> buf;    y = std::stod(buf, &sz);
        file >> buf;    perm = visc * exp(std::stod(buf, &sz));
        file >> buf;
        id_x = 1 + (x - x1) / hx;
        id_y = 1 + (y - y1) / hy;
        //if (id != 17)
        //{
            wells.push_back(Well(id, (num_y + 2) * id_x + id_y));
            conds.push_back({ (num_y + 2) * id_x + id_y, perm });
        //}
    }

    file.close();
}
int main()
{
    double x1 = 8402.8, x2 = 13900.0;// x2 = 13988.4;
    double y1 = 24917.4, y2 = 29700.0;// y2 = 30242.5;
    int num_x = 41, num_y = 41;

	dual_stoch_oil::Properties props;
	
	props.possible_steps_num = 4;
	props.start_time_simple_approx = 2;
	props.t_dim = 3600.0;
	props.ht = props.ht_min = 10000000.0;
	props.ht_max = 100000000.0;
    props.hx = props.R_dim = x2 - x1;
    props.hy = y2 - y1;
    props.hz = 10.0;
	props.num_x = num_x;
    props.num_y = num_y;
	int num = (props.num_x + 2) * (props.num_y + 2);

	props.props_sk.p_init = props.props_sk.p_out = 275.39 * BAR_TO_PA;
	props.props_sk.perm = 500.0;
	props.props_sk.m = 0.1;
	props.props_sk.beta = 4.E-10;
    props.props_sk.l_f = 500.0;
	props.props_sk.sigma_f = 0.66;

	props.props_oil.visc = 1.0;
	props.props_oil.rho_stc = 887.261;
	props.props_oil.beta = 1.0 * 1.e-9;
	props.props_oil.p_ref = props.props_sk.p_init;

	//props.wells.push_back(Well(0, (props.num_y + 2) * (int)(props.num_x / 2 + 1) + (int)(props.num_x / 2 + 1)));
    loadWells(x1, x2, y1, y2, num_x, num_y, "props/wells_gen.txt", props.wells, props.conditions, props.props_oil.visc);
    for (auto& well1 : props.wells)
    {
        well1.periodsNum = 1;
        well1.period.resize(well1.periodsNum);
        well1.period[0] = 365.0 * 86400.0;
        well1.rate.resize(well1.periodsNum);
        //if(well1.id != 13)
            well1.rate[0] = -1000.0;
        //else
        //    well1.rate[0] = -100.0;
        well1.pwf.resize(well1.periodsNum);
        //well1.pwf[0] = 180.0 * BAR_TO_PA;
        well1.leftBoundIsRate.resize(well1.periodsNum);
        well1.leftBoundIsRate[0] = true;
        well1.rw = 0.1;
    }

	Scene<issues::DualStochOil> scene;
	scene.load(props);
	scene.start();

	return 0;
}