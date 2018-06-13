#ifndef SCENE_HPP_
#define SCENE_HPP_

#include <memory>
#include "paralution.hpp"

template <typename issueType>
class Scene
{
public:
	typename typedef issueType::Model Model;
	typename typedef Model::Properties Properties;
	typename typedef issueType::Method Method;
protected:
	std::shared_ptr<Model> model;
	std::shared_ptr<Method> method;
public:

	Scene() {};
	~Scene() { paralution::stop_paralution(); };

	void load(const Properties& props)
	{
		paralution::init_paralution();

		model = std::make_shared<Model>();
		model->load(props);
		model->setSnapshotter(model.get());
		method = std::make_shared<Method>(model.get());

		model->snapshot_all(0);
	}
	void start()
	{
		method->start();
	};
};

#endif /* SCENE_HPP_ */
