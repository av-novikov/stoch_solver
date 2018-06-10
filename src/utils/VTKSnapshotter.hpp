#ifndef VTKSNAPSHOTTER_HPP_
#define VTKSNAPSHOTTER_HPP_

namespace snapshotter
{
	template<class modelType>
	class VTKSnapshotter
	{
	public:
		typedef modelType Model;
		typedef typename Model::Mesh Mesh;
		typedef typename Model::Cell Cell;
	protected:
		const Model* model;
		const Mesh* mesh;
		const std::string prefix = "snaps/";
		std::string pattern;
		std::string replace(std::string filename, std::string from, std::string to);
		std::string getFileName(const int snap_idx);

		double R_dim;
		size_t num_x, num_y;
	public:
		VTKSnapshotter(const Model* _model);
		~VTKSnapshotter();

		void dump(const int i);
	};
};

#endif /* VTKSNAPSHOTTER_HPP_ */