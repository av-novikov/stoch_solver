#ifndef ABSTRACTMODEL_HPP_
#define ABSTRACTMODEL_HPP_

#include <vector>
#include <string>
#include <map>
#include <memory>

#include "src/utils/VTKSnapshotter.hpp"
#include "src/Well.hpp"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

template <typename TVarContainer, typename propsType, template <typename TVarContainer> class TVariables, class MeshType, class modelType>
class AbstractModel : public TVariables<TVarContainer>
{
	template<typename> friend class snapshotter::VTKSnapshotter;
	template<typename> friend class AbstractSolver;
public:
	typedef TVarContainer VarContainer;
	typedef TVariables<TVarContainer> Variables;
	typedef MeshType Mesh;
	typename typedef Mesh::Cell Cell;
	typedef propsType Properties;
	typedef typename snapshotter::VTKSnapshotter<modelType> SnapShotter;
protected:
	std::shared_ptr<Mesh> mesh;
	std::shared_ptr<SnapShotter> snapshotter;

	// Spacial properties
	double r_w;
	double Volume;
	size_t cellsNum;
	// Number of unknown variables
	size_t varNum;
	// Number of unknown variables in cell
	static const int var_size;

	// Rate of the well
	std::vector<Well> wells;

	// Temporary properties
	double ht;
	double ht_min;
	double ht_max;
	// Number of periods
	size_t periodsNum;
	// End times of periods [sec]
	std::vector<double> period;
	// Oil rates [m3/day]
	std::vector<double> rate;
	// Vector of BHPs [bar]
	std::vector<double> pwf;
	// If left boundary condition would be 2nd type
	bool leftBoundIsRate;
	// If right boundary condition would be 1st type
	bool rightBoundIsPres;
	// BHP will be converted to the depth
	double depth_point;
	// During the time flow rate decreases 'e' times in well test [sec] 

	virtual void loadMesh(const std::string nebrFileName)
	{
		/*mshreader::MshReader reader;
		mesh = std::make_shared<grid::Mesh>(*reader.read(nebrFileName, R_dim));
		mesh->process_geometry();

		cellsNum = mesh.get()->getCalcCellsSize();
		varNum = VarContainer::size * cellsNum;
		u_prev.resize(varNum);
		u_iter.resize(varNum);
		u_next.resize(varNum);
		Volume = mesh->Volume;*/
	}
	adouble linearInterp1d(const adouble a1, const double r1, const adouble a2, const double r2)
	{
		return (a1 * r2 + a2 * r1) / (r1 + r2);
	};
	virtual void setProps(const propsType& props) = 0;
	virtual void makeDimLess() = 0;
	virtual void setInitialState() = 0;
public:
	AbstractModel()
	{
		grav = 9.8;
	};
	virtual ~AbstractModel() {};

	// Dimensions
	double t_dim;
	double R_dim;
	double Z_dim;
	double P_dim;
	double T_dim;
	double Q_dim;
	double grav;

	void load(const Properties& props, const std::string nebrFileName)
	{
		setProps(props);
		loadMesh(nebrFileName);
		setInitialState();
	};
	virtual void setPeriod(const int period) = 0;
	virtual void setWellborePeriod(int period, double cur_t) {};
	int getCellsNum() { return cellsNum; };
	void snapshot_all(const int i) { snapshotter->dump(i); }
	const Mesh* getMesh() const
	{
		return mesh.get();
	}
	Mesh* getMesh()
	{
		return mesh.get();
	}
	void setSnapshotter(const modelType* model)
	{
		snapshotter = std::make_shared<SnapShotter>(model);
	}
};

#endif /* ABSTRACTMODEL_HPP_ */
