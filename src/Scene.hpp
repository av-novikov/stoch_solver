#ifndef SCENE_HPP_
#define SCENE_HPP_

#include <memory>

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
	~Scene() {};

	void load(const Properties& props)
	{
		model = std::make_shared<Model>();
		model->setProps(props);
		model->setSnapshotter(model.get());
		//method = std::make_shared<Method>(model.get());
		method = std::make_shared<Method>();

		model->snapshot_all(0);
	}
};

#endif /* SCENE_HPP_ */
