#ifndef DUAL_STOCH_OIL_HPP_
#define DUAL_STOCH_OIL_HPP_

#include "src/model/AbstractModel.hpp"
#include "src/grid/Variables.hpp"
#include "src/grid/Mesh.hpp"
#include "src/model/dual_stoch_oil/Properties.hpp"
#include "src/Well.hpp"
#include "paralution.hpp"

namespace dual_stoch_oil
{
	/*typedef var::containers::TapeVar1Phase TapeVariable0;
	typedef var::containers::TapeStochVar1Phase1 TapeVariable1;
	typedef var::containers::TapeStochVar1Phase2 TapeVariable2;
	typedef var::containers::TapeStochVar1Phase3 TapeVariable3;*/
	class DualStochOil : public AbstractDualGridModel<Properties,mesh::CellRectangularUniformGrid, mesh::NodeRectangularUniformGrid,DualStochOil,var::DualStochVariables,
				var::containers::Var1phase>
	{
		template<typename> friend class snapshotter::VTKSnapshotter;
		template<typename> friend class AbstractDualGridMethod;
		friend class DualStochOilMethod;
	public:

	protected:
		void makeDimLess();
		void setInitialState();

		adouble *x_cell, *x_node;
		adouble *h_cell, *h_node;

		int possible_steps_num, start_time_simple_approx;
		Skeleton_Props props_sk;
		Oil_Props props_oil;
		std::vector<Well> wells;
        // Conditioning = Kriging
        std::vector<Measurement> conditions;
        double* inv_cond_cov;
        std::vector<double> Favg_cells, Favg_nodes;
        std::vector<std::vector<double>> Cf_cells, Cf_nodes;

        void loadPermAvg(const std::string fileName);
        void writeCPS(const int i);
        template<class TElem>
		inline double getPerm_prior(const TElem& elem) const
		{
            /*int ind_x = int(cell.id / (mesh->num_y + 2));
            int ind_y = cell.id % (mesh->num_y + 2);
            int id = cell.id;
            if (cell.type != elem::QUAD)
            {
                //if ((ind_x == 0 || ind_x == mesh->num_x + 1) && (ind_y == 0 || ind_y == mesh->num_y + 1))
                //    return 1.0;
                //else
                //{
                    if (ind_y == 0)
                        id += 1;
                    else if (ind_y == mesh->num_y + 1)
                        id -= 1;
                    if (ind_x == 0)
                        id += mesh->num_y + 2;
                    else if (ind_x == mesh->num_x + 1)
                        id -= mesh->num_y + 2;
                //}
            }
            ind_x = int(id / (mesh->num_y + 2));
            ind_y = mesh->num_y + 1 - id % (mesh->num_y + 2);
            return props_sk.perm_grd[(ind_x - 1) * mesh->num_y + (ind_y - 1)];*/
            return props_sk.perm;
		};
        template<class TElem>
        inline double getS(const TElem& cell) const
        {
            return props_sk.m * (props_oil.beta + props_sk.beta);
        };
        // Apriori (unconditioned) statistical moments
        template<class TElem>
        inline double getFavg_prior(const TElem& elem) const
        {
            return log(getPerm_prior(elem) / props_oil.visc) - getSigma2f_prior(elem) / 2.0;
        };
        inline double getCf_coord(const point::Point& p1, const point::Point& p2) const
        {
            const double dist = point::distance(p1, p2) / props_sk.l_f;
            return props_sk.sigma_f * props_sk.sigma_f * exp(-dist * dist);
            //return props_sk.sigma_f * props_sk.sigma_f * exp(-dist);
        };
        template<class TElem>
		inline double getCf_prior(const TElem& elem1, const TElem& elem2) const
		{
			return getCf_coord(elem1.cent, elem2.cent);
		};
        template<class TElem>
		inline double getSigma2f_prior(const TElem& elem) const
		{
			return getCf_prior(elem, elem);
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
                    const auto& cell1 = cell_mesh->cells[conditions[i].id];
                    for (int j = 0; j < conditions.size(); j++)
                    {
                        const auto& cell2 = cell_mesh->cells[conditions[j].id];

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
                /*for (int i = 0; i < conditions.size(); i++)
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
                // Mult for cell grid
                std::vector<std::vector<double>> cell_mult_mat;
                cell_mult_mat.resize(cellsNum);
                std::for_each(cell_mult_mat.begin(), cell_mult_mat.end(), [&](std::vector<double>& vec) { vec.resize(conditions.size(), 0.0); });
                for (int i = 0; i < cellsNum; i++)
                {
                    const Cell& cell = cell_mesh->cells[i];
                    for (int k = 0; k < conditions.size(); k++)
                    {
                        const Cell& c_cell = cell_mesh->cells[conditions[k].id];
                        s = 0.0;
                        for (int k1 = 0; k1 < conditions.size(); k1++)
                        {
                            const Cell& c_cell1 = cell_mesh->cells[conditions[k1].id];
                            s += getCf_prior(cell, c_cell1) * inv_cond_cov[k1 * conditions.size() + k];
                        }
                        cell_mult_mat[i][k] = s;
                    }
                }
                // Mult for node grid
                std::vector<std::vector<double>> node_mult_mat;
                node_mult_mat.resize(nodesNum);
                std::for_each(node_mult_mat.begin(), node_mult_mat.end(), [&](std::vector<double>& vec) { vec.resize(conditions.size(), 0.0); });
                for (int i = 0; i < nodesNum; i++)
                {
                    const Node& node = node_mesh->nodes[i];
                    for (int k = 0; k < conditions.size(); k++)
                    {
                        const Node& c_node = node_mesh->nodes[conditions[k].id];
                        s = 0.0;
                        for (int k1 = 0; k1 < conditions.size(); k1++)
                        {
                            const Node& c_node1 = node_mesh->nodes[conditions[k1].id];
                            s += getCf_prior(node, c_node1) * inv_cond_cov[k1 * conditions.size() + k];
                        }
                        node_mult_mat[i][k] = s;
                    }
                }
                // Kriging for cells
                for (int i = 0; i < cellsNum; i++)
                {
                    const Cell& cell1 = cell_mesh->cells[i];
                    for (int k2 = 0; k2 < conditions.size(); k2++)
                    {
                        const auto& cond = conditions[k2];
                        const auto c_cell = cell_mesh->cells[cond.id];
                        Favg_cells[i] += cell_mult_mat[i][k2] * (log(cond.perm / props_oil.visc) - getFavg_prior(c_cell));
                    }
                    // Covariance
                    for (int j = 0; j < cellsNum; j++)
                    {
                        const Cell& cell2 = cell_mesh->cells[j];
                        for (int k2 = 0; k2 < conditions.size(); k2++)
                        {
                            const auto& cond = conditions[k2];
                            Cf_cells[i][j] -= cell_mult_mat[i][k2] * getCf_prior(cell_mesh->cells[cond.id], cell2);
                        }
                    }
                    if (Cf_cells[i][i] < 0.0 && Cf_cells[i][i] > -EQUALITY_TOLERANCE)
                        Cf_cells[i][i] = 0.0;
                }
                // Kriging for nodes
                for (int i = 0; i < nodesNum; i++)
                {
                    const Node& node1 = node_mesh->nodes[i];
                    for (int k2 = 0; k2 < conditions.size(); k2++)
                    {
                        const auto& cond = conditions[k2];
                        const auto c_node = node_mesh->nodes[cond.id];
                        Favg_nodes[i] += node_mult_mat[i][k2] * (log(cond.perm / props_oil.visc) - getFavg_prior(c_node));
                    }
                    // Covariance
                    for (int j = 0; j < nodesNum; j++)
                    {
                        const Node& node2 = node_mesh->nodes[j];
                        for (int k2 = 0; k2 < conditions.size(); k2++)
                        {
                            const auto& cond = conditions[k2];
                            Cf_nodes[i][j] -= node_mult_mat[i][k2] * getCf_prior(node_mesh->nodes[cond.id], node2);
                        }
                    }
                    if (Cf_nodes[i][i] < 0.0 && Cf_nodes[i][i] > -EQUALITY_TOLERANCE)
                        Cf_nodes[i][i] = 0.0;
                }
            }

            for (const auto& cond : conditions)
            {
                auto it = find_if(wells.begin(), wells.end(), [&](const Well& well) {return well.cell_id == cond.id; });
                if (it != wells.end())
                    it->isCond = true;
            }

            for (auto& well : wells)
            {
                const auto& cell = cell_mesh->cells[well.cell_id];
                well.perm = exp(getFavg(cell) + getSigma2f(cell) / 2.0) * props_oil.visc;
            }
        };
        template<class TElem>
        inline double getPerm(const TElem& elem) const
        {
            return props_oil.visc * exp(getFavg(elem) + getSigma2f(elem) / 2.0);
        };
        template<class TElem>
        inline double getFavg(const TElem& elem) const
        {
            return 0.0;
        };
        inline double getFavg(const Cell& elem) const
        {
            return Favg_cells[elem.id];
        };
        inline double getFavg(const Node& elem) const
        {
            return Favg_nodes[elem.id];
        };
        template<class TElem>
        inline double getCf(const TElem& elem1, const TElem& elem2) const
        {
            return 0.0;
        };
        inline double getCf(const Cell& elem1, const Cell& elem2) const
        {
            return Cf_cells[elem1.id][elem2.id];
        };
        inline double getCf(const Node& elem1, const Node& elem2) const
        {
            return Cf_nodes[elem1.id][elem2.id];
        };
        template<class TElem>
        inline double getSigma2f(const TElem& elem) const
        {
            return getCf(elem, elem);
        };
        template<class TElem>
        inline double getKg(const TElem& elem) const
        {
            return exp(getFavg(elem));
        };
        template<class TElem>
        inline double getGeomPerm(const TElem& elem) const
        {
            return getKg(elem) * props_oil.visc;
        };

		adouble solveInner_p0(const Cell& cell) const;
		adouble solveBorder_p0(const Cell& cell) const;
		adouble solveSource_p0(const Well& well) const;

		adouble solveInner_Cfp(const Cell& cell, const Cell& cur_cell) const;
		adouble solveBorder_Cfp(const Cell& cell, const Cell& cur_cell) const;
		adouble solveSource_Cfp(const Well& well, const Cell& cur_cell) const;

        adouble solveInnerNode_Cfp(const Node& node, const Node& cur_node) const;
        adouble solveBorderNode_Cfp(const Node& node, const Node& cur_node) const;

		adouble solveInner_p2(const Cell& cell) const;
		adouble solveBorder_p2(const Cell& cell) const;
		adouble solveSource_p2(const Well& well) const;

		adouble solveInner_Cp(const Cell& cell, const Cell& cur_cell, const size_t step_idx, const size_t cur_step_idx) const;
		adouble solveBorder_Cp(const Cell& cell, const Cell& cur_cell, const size_t step_idx) const;
		adouble solveSource_Cp(const Well& well, const Cell& cur_cell, const size_t step_idx) const;

        double getRate(const Well& well) const;
        double getRateVar(const Well& well, const int step_idx) const;
        double getPwf(const Well& well) const;
        double getPwfVar(const Well& well, const int step_idx) const;
	public:
		DualStochOil();
		~DualStochOil();

		void setProps(const Properties& props);
		void setPeriod(const int period);
	};
};

#endif /* DUAL_STOCH_OIL_HPP_ */