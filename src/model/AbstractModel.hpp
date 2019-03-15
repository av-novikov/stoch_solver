#ifndef ABSTRACTMODEL_HPP_
#define ABSTRACTMODEL_HPP_

#include <vector>
#include <string>
#include <map>
#include <memory>

#include "src/utils/VTKSnapshotter.hpp"
#include "src/Well.hpp"
#include "src/grid/Variables.hpp"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

template <typename propsType,
			class MeshType,
			class modelType,
			template <typename TVarContainer> class TVariables,
			typename TVarContainer>
class AbstractModel : public TVariables<TVarContainer>
{
	template<typename> friend class snapshotter::VTKSnapshotter;
	template<typename> friend class AbstractMethod;
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
	static const int var_size = VarContainer::size;

	// Rate of the well
	std::vector<Well> wells;

	// Temporary properties
	double ht;
	double ht_min;
	double ht_max;
	// If right boundary condition would be 1st type
	bool rightBoundIsPres;
	// BHP will be converted to the depth
	double depth_point;
	// During the time flow rate decreases 'e' times in well test [sec] 

	adouble linearInterp1d(const adouble a1, const double r1, const adouble a2, const double r2) const
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

	void load(const Properties& props)
	{
		setProps(props);
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
	virtual double getRate(const Well& well) const = 0;
};

template <typename propsType,
    class MeshType1, class MeshType2,
    class modelType,
    template <typename TVarContainer> class TVariables,
    typename TVarContainer>
    class AbstractDualGridModel : public TVariables<TVarContainer>
{
    template<typename> friend class snapshotter::VTKSnapshotter;
    template<typename> friend class AbstractMethod;
public:
    typedef TVarContainer VarContainer;
    typedef TVariables<TVarContainer> Variables;
    typedef MeshType1 CellMesh;
    typedef MeshType2 NodeMesh;
    typedef CellMesh Mesh;
    typename typedef CellMesh::Cell Cell;
    typename typedef NodeMesh::Node Node;
    typedef propsType Properties;
    typedef typename snapshotter::VTKSnapshotter<modelType> SnapShotter;
protected:
    std::shared_ptr<CellMesh> cell_mesh;
    std::shared_ptr<NodeMesh> node_mesh;
    std::shared_ptr<SnapShotter> snapshotter;

    // Spacial properties
    double r_w;
    double Volume;
    int cellsNum, nodesNum;
    // Number of unknown variables
    size_t varNum;
    // Number of unknown variables in cell
    static const int var_size = VarContainer::size;

    // Rate of the well
    std::vector<Well> wells;
     
    // Temporary properties
    double ht;
    double ht_min;
    double ht_max;
    // If right boundary condition would be 1st type
    bool rightBoundIsPres;
    // BHP will be converted to the depth
    double depth_point;
    // During the time flow rate decreases 'e' times in well test [sec] 

    adouble linearInterp1d(const adouble a1, const double r1, const adouble a2, const double r2) const
    {
        return (a1 * r2 + a2 * r1) / (r1 + r2);
    };
    virtual void setProps(const propsType& props) = 0;
    virtual void makeDimLess() = 0;
    virtual void setInitialState() = 0;
public:
    AbstractDualGridModel()
    {
        grav = 9.8;
    };
    virtual ~AbstractDualGridModel() {};

    // Dimensions
    double t_dim;
    double R_dim;
    double Z_dim;
    double P_dim;
    double T_dim;
    double Q_dim;
    double grav;

    void load(const Properties& props)
    {
        setProps(props);
        setInitialState();
    };
    virtual void setPeriod(const int period) = 0;
    virtual void setWellborePeriod(int period, double cur_t) {};
    int getCellsNum() { return cellsNum; };
    void snapshot_all(const int i) { snapshotter->dump(i); }
    const CellMesh* getCellMesh() const
    {
        return cell_mesh.get();
    }
    CellMesh* getCellMesh()
    {
        return cell_mesh.get();
    }
    const NodeMesh* getNodeMesh() const
    {
        return node_mesh.get();
    }
    NodeMesh* getNodeMesh()
    {
        return node_mesh.get();
    }
    void setSnapshotter(const modelType* model)
    {
        snapshotter = std::make_shared<SnapShotter>(model);
    }
    virtual double getRate(const Well& well) const = 0;
};

#endif /* ABSTRACTMODEL_HPP_ */
