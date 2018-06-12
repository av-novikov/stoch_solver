#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/oil/Properties.hpp"
#include "src/Well.hpp"

namespace oil
{
	typedef var::containers::TapeVar1Phase TapeVariable;
	class Oil : public AbstractModel<var::containers::Var1phase, Properties, var::BasicVariables, mesh::RectangularUniformGrid, Oil>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class OilMethod;
	protected:
		void makeDimLess();
		void setInitialState();

		TapeVariable* x;
		adouble* h;

		Skeleton_Props props_sk;
		Oil_Props props_oil;
		std::vector<Well> wells;
	public:
		Oil();
		~Oil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
	};
};