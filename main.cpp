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
	
	props.hx = props.hy = props.R_dim = 100.0;		props.hz = 1.0;
	props.num_x = props.num_y = 10;

	Scene<issues::Oil> scene;
	scene.load(props);


	return 0;
}