#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/oil/Properties.hpp"

namespace oil
{
	typedef var::containers::TapeVar1Phase TapeVariable;
	class Oil : public AbstractModel<var::containers::Var1phase, Properties, var::BasicVariables, mesh::RectangularUniformGrid, Oil>
	{
	protected:
		void makeDimLess();
		void setInitialState();
	public:
		Oil();
		~Oil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
	};
};