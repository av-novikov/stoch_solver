#ifndef STOCH_OIL_HPP_
#define STOCH_OIL_HPP_

#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/RectangularUniformGrid.hpp"
#include "src/model/stoch_oil/Properties.hpp"
#include "src/Well.hpp"
#include "paralution.hpp"

namespace stoch_oil
{
	/*typedef var::containers::TapeVar1Phase TapeVariable0;
	typedef var::containers::TapeStochVar1Phase1 TapeVariable1;
	typedef var::containers::TapeStochVar1Phase2 TapeVariable2;
	typedef var::containers::TapeStochVar1Phase3 TapeVariable3;*/
	class StochOil : public AbstractModel<Properties,mesh::RectangularUniformGrid,StochOil,var::StochVariables,
				var::containers::Var1phase>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class StochOilMethod;
	public:

	protected:
		void makeDimLess();
		void setInitialState();

		adouble* x;
		adouble* h;

		int possible_steps_num, start_time_simple_approx;
		Skeleton_Props props_sk;
		Oil_Props props_oil;
		std::vector<Well> wells;
        // Conditioning = Kriging
        std::vector<Measurement> conditions;
        double* inv_cond_cov;
        std::vector<double> Favg;
        std::vector<std::vector<double>> Cf;

        void loadPermAvg(const std::string fileName);
        void writeCPS(const int i);
		inline double getPoro(const Cell& cell) const
		{
			return props_sk.m;
		};
		inline double getPerm(const Cell& cell) const
		{
			return props_sk.perm;
            //const int id = (cell.type == elem::QUAD) ? cell.id : cell.stencil[1];
            //return props_sk.perm_grd[(id / (mesh->num_y + 2) - 1) * mesh->num_y + id % (mesh->num_y + 2) - 1];
		};
        inline double getS(const Cell& cell) const
        {
            return getPoro(cell) * (props_oil.beta + props_sk.beta);
        };

        // Apriori (unconditioned) statistical moments
        inline double getFavg_prior(const Cell& cell) const
        {
            return log(getPerm(cell) / props_oil.visc);
        };
		inline double getCf_prior(const Cell& cell, const Cell& beta) const
		{
			auto corrFoo1 = [this](const auto& p1, const auto& p2) -> double
			{
				return props_sk.sigma_f * props_sk.sigma_f * exp(-point::distance(p1, p2) / props_sk.l_f);
			};
			auto corrFoo2 = [this](const auto& p1, const auto& p2) -> double
			{
				const double dist = point::distance(p1, p2) / props_sk.l_f;
				return props_sk.sigma_f * props_sk.sigma_f * exp(-dist * dist);
			};

			return corrFoo1(cell.cent, beta.cent);
		};
		inline double getSigma2f_prior(const Cell& cell) const
		{
			return getCf(cell, cell);
		};
        // Apostreriori (conditioned) statistical moments
        void calculateConditioning()
        {
            const int matSize = conditions.size() * conditions.size();

            if (matSize > 0)
            {
                int* ind_i = new int[matSize];
                int* ind_j = new int[matSize];
                double* cond_cov = new double[matSize];

                int counter = 0;
                for (int i = 0; i < conditions.size(); i++)
                {
                    const auto& cell1 = mesh->cells[conditions[i].id];
                    for (int j = 0; j < conditions.size(); j++)
                    {
                        const auto& cell2 = mesh->cells[conditions[j].id];

                        ind_i[counter] = i;     ind_j[counter] = j;
                        cond_cov[counter] = getCf_prior(cell1, cell2);
                        counter++;
                    }
                }

                paralution::LocalMatrix<double> A;
                paralution::Inversion<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double> ls;
                A.Assemble(ind_i, ind_j, cond_cov, matSize, "A", conditions.size(), conditions.size());
                ls.SetOperator(A);
                ls.Build();
                ls.inverse_.LeaveDataPtrDENSE(&inv_cond_cov);
                ls.Clear();

                double s;
                // Check inversion
                /*counter = 0;
                for (int i = 0; i < conditions.size(); i++)
                {
                    for (int j = 0; j < conditions.size(); j++)
                    {
                        s = 0.0;
                        for (int k = 0; k < conditions.size(); k++)
                        {
                            double q1 = cond_cov[i * conditions.size() + k];
                            double q2 = inv_cond_cov[k * conditions.size() + j];
                             s += q1 * q2;
                        }
                        if(i == j)
                            assert(fabs(s - 1) < EQUALITY_TOLERANCE);
                        else
                            assert(fabs(s) < EQUALITY_TOLERANCE);
                    }
                }*/
                delete[] ind_i, ind_j;
                delete[] cond_cov;

                std::vector<std::vector<double>> mult_mat;
                mult_mat.resize(cellsNum);
                std::for_each(mult_mat.begin(), mult_mat.end(), [&](std::vector<double>& vec) { vec.resize(conditions.size(), 0.0); });
                for (int i = 0; i < cellsNum; i++)
                {
                    const Cell& cell = mesh->cells[i];
                    for (int k = 0; k < conditions.size(); k++)
                    {
                        const Cell& c_cell = mesh->cells[conditions[k].id];
                        s = 0.0;
                        for (int k1 = 0; k1 < conditions.size(); k1++)
                        {
                            const Cell& c_cell1 = mesh->cells[conditions[k1].id];
                            s += getCf_prior(cell, c_cell1) * inv_cond_cov[k1 * conditions.size() + k];
                        }
                        mult_mat[i][k] = s;
                    }
                }

                for (int i = 0; i < cellsNum; i++)
                {
                    const Cell& cell1 = mesh->cells[i];
                    for (int k2 = 0; k2 < conditions.size(); k2++)
                    {
                        const auto& cond = conditions[k2];
                        Favg[i] += mult_mat[i][k2] * (log(cond.perm / props_oil.visc) - getFavg_prior(mesh->cells[cond.id]));
                    }
                    // Covariance
                    for (int j = 0; j < cellsNum; j++)
                    {
                        const Cell& cell2 = mesh->cells[j];
                        for (int k2 = 0; k2 < conditions.size(); k2++)
                        {
                            const auto& cond = conditions[k2];
                            Cf[i][j] -= mult_mat[i][k2] * getCf_prior(mesh->cells[cond.id], cell2);
                        }
                    }
                }
            }

            for (auto& well : wells)
                well.perm = exp(getFavg(mesh->cells[well.cell_id])) * props_oil.visc;
        };
        inline double getFavg(const Cell& cell) const
        {
            return Favg[cell.id];
        };
        inline double getCf(const Cell& cell, const Cell& beta) const
        {
            return Cf[cell.id][beta.id];
        };
        inline double getSigma2f(const Cell& cell) const
        {
            return getCf(cell, cell);
        };
        inline double getKg(const Cell& cell) const
        {
            return exp(getFavg(cell));
        };

		adouble solveInner_p0(const Cell& cell) const;
		adouble solveBorder_p0(const Cell& cell) const;
		adouble solveSource_p0(const Well& well) const;

		adouble solveInner_Cfp(const Cell& cell, const Cell& cur_cell) const;
		adouble solveBorder_Cfp(const Cell& cell, const Cell& cur_cell) const;
		adouble solveSource_Cfp(const Well& well, const Cell& cur_cell) const;

		adouble solveInner_p2(const Cell& cell) const;
		adouble solveBorder_p2(const Cell& cell) const;
		adouble solveSource_p2(const Well& well) const;

		adouble solveInner_Cp(const Cell& cell, const Cell& cur_cell, const size_t step_idx, const size_t cur_step_idx) const;
		adouble solveBorder_Cp(const Cell& cell, const Cell& cur_cell, const size_t step_idx) const;
		adouble solveSource_Cp(const Well& well, const Cell& cur_cell, const size_t step_idx) const;
	public:
		StochOil();
		~StochOil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
		double getRate(const Well& well) const;
		double getPwf(const Well& well) const;
	};
};

#endif /* STOCH_OIL_HPP_ */