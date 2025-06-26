#pragma once

#include "util.hpp"
#include "Component.hpp"

class Node
{
    public:
        static float cvx_thresh_;
        static std::function<float(const Node*)> m_g;
        static std::function<float(const Node*)> m_f;
        Node(
            const std::vector<Component*> &comps, 
            const LCC_3 &lcc,
            const int max_possible_comps,
            const float initial_volume_ch,
            const size_t depth,
            const float sigma,
            const std::pair<float, float> &sym_axis = {-1, 0}): 
            m_comps(comps), 
            m_max_possible_comps(max_possible_comps),
            m_initial_volume_ch(initial_volume_ch),
            m_depth(depth),
            m_sigma(sigma),
            m_sym_axis(sym_axis)
            {
                m_lcc = lcc;
                update_volume();
                update_convexity();
                update_degree();
                if (sym_axis.first >= 0)
                    update_IoS();
            };
        
        
        ~Node()
        {
            for (auto &c : m_comps)
            {
                delete c;
                c = nullptr;
            }
        };

        float g() const {return m_g(this);};
        float h() const {return m_f(this);};
        size_t depth() const {return m_depth;};

        int num_potential_cuts() const 
        {
            int num_cuts = 0;
            for (auto comp : m_comps)
                num_cuts += comp->num_potential_cuts();
            return num_cuts;
        }

        void update_volume()
        {
            m_sum_volume = 0;
            m_volumes.clear();
            m_volumes.resize(m_comps.size());
            m_sum_delta_ch = 0;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                m_volumes[i] = m_comps[i]->volume();
                m_sum_volume += m_volumes[i];
                m_sum_delta_ch += m_comps[i]->delta_ch();
            }
            m_initial_delta_ch = m_initial_volume_ch - m_sum_volume;
        }

        void update_convexity()
        {
            m_is_convex.clear();
            m_is_convex.resize(m_comps.size());
            for (size_t i = 0; i < m_comps.size(); i++)
                m_is_convex[i] = m_comps[i]->num_cc_edges() == 0;
        }

        int num_comps() const {return m_comps.size();};

        #if VIS_NODE == 1 || VIS_CHILD == 1 || VIS == 1
        void draw(std::string title = "") const {
            std::string full_tilte = title + " ";
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                title += std::to_string(m_comps[i]->num_cc_edges()) 
                        + "(" + std::to_string(m_comps[i]->num_potential_cuts()) + ")" +"-";
            }
            title += "(d="+ std::to_string(this->reduced_delta_ch()) 
                     + "=" + std::to_string(m_initial_delta_ch) + "-" + std::to_string(m_sum_delta_ch);
            title += ", rcc=" + std::to_string(remaining_ccv_edges()) + "/" + std::to_string(m_max_possible_comps);
            title += " c=" + std::to_string(this->cost()) + ", h=" + std::to_string(this->h()) 
                  + ", g=" + std::to_string(this->g()) + ", d=" + std::to_string(m_depth) + ")";
            CGAL::draw(m_lcc, title.c_str());
            
        };
        #endif

        bool operator== (const Node &other) const
        {
            std::vector<float> other_volumes = other.m_volumes;
            std::vector<float> this_volumes = this->m_volumes;
            if (this_volumes.size() != other_volumes.size())
                return false;
            std::sort(this_volumes.begin(), this_volumes.end());
            std::sort(other_volumes.begin(), other_volumes.end());
            for (size_t i = 0; i < this_volumes.size(); i++)
            {
                if (std::abs(this_volumes[i] - other_volumes[i]) > Global_var::same_vol_factor * m_sum_volume)
                    return false;
            }
            // #if VIS==1
            // this->draw("this");
            // other.draw("other");
            // #endif
            return true;
        };

        bool operator< (const Node &other) const
        {
            if (*this == other) return false;
            
            if (std::abs(this->cost() - other.cost()) > 1e-3)
                return this->cost() < other.cost();
            else
            {
                // if (std::abs(this->h() - other.h()) > 1e-3) return this->h() < other.h();
                // else 
                // {
                //     return true; // preserve the node
                // }
                return this->h() < other.h();
            }
        };

        float cost() const {
            return g() + m_sigma * h();
        };

        void get_comp_cut_combinations(
            Eigen::MatrixXi &combinations) const
        {
            int num_all_combs = 1;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (m_is_convex[i] || m_comps[i]->num_potential_cuts() == 0) continue;
                num_all_combs *= m_comps[i]->num_potential_cuts();
            }
            combinations = Eigen::MatrixXi::Constant(num_all_combs, m_comps.size(), -1);

            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (m_is_convex[i] || m_comps[i]->num_potential_cuts() == 0) continue;
                int remaining_num_combs = 1;
                for (size_t j = i + 1; j < m_comps.size(); j++)
                {
                    if (m_is_convex[j] || m_comps[j]->num_potential_cuts() == 0) continue;
                    remaining_num_combs *= m_comps[j]->num_potential_cuts();
                }
                int num_sub_combinations = remaining_num_combs * m_comps[i]->num_potential_cuts();
                int num_round = num_all_combs / num_sub_combinations;
                for (int k = 0; k < num_round; k++)
                    for (int j = 0; j < m_comps[i]->num_potential_cuts(); j++)
                        combinations.block(k * num_sub_combinations + j * remaining_num_combs, i, remaining_num_combs, 1) 
                                    = Eigen::MatrixXi::Constant(remaining_num_combs, 1, j);
            }
        }

        float initial_volume_ch() const
        { return m_initial_volume_ch;}

        void expand(
            std::vector<std::shared_ptr<const Node>> &children) const 
        {
            assert (m_is_convex.size() == m_comps.size());
            // std::vector<size_t> comp_to_cut = sort_comps_to_cut();
            // if (comp_to_cut.size() == 0)
            // {
            //     std::cout << "no comp to cut\n";
            //     return;
            // }
            
            // parallelize this loop
            // FAST EXPAND
            Eigen::MatrixXi combinations;
            get_comp_cut_combinations(combinations);

            std::vector<std::shared_ptr<const Node>> node_children;
            node_children.resize(combinations.rows());
            #pragma omp parallel for
            for (size_t i = 0; i < combinations.rows(); i++)
            {
                LCC_3 lcc_copy = m_lcc;
                lcc_copy.onsplit_functor<1>()=Split_functor_1(lcc_copy);
                lcc_copy.onsplit_functor<2>()=Split_functor_2(lcc_copy);
                lcc_copy.onsplit_functor<3>()=Split_functor_3(lcc_copy);

                bool new_cut = false;
                std::vector<Component*> new_comps;
                for (size_t j = 0; j < m_comps.size(); j++)
                {
                    if (combinations(i, j) < 0)
                        new_comps.push_back(new Component(*m_comps[j]));
                    else
                    {
                        Component *comp = m_comps[j];
                        const Plane* cut = comp->get_cut(combinations(i, j));
                        std::vector<Component*> cut_comps;
                        std::cout << "cutting comp " << j << " by cut " << combinations(i, j) << '\n';
                        int result = comp->cut(cut, lcc_copy, cut_comps);
                        if (result == EXIT_SUCCESS)
                        {
                            new_comps.insert(new_comps.end(), cut_comps.begin(), cut_comps.end());
                            new_cut = true;
                        }
                        else
                            new_comps.push_back(new Component(*comp));
                    }
                }
                if (new_cut)
                {
                    Node* new_node = new Node(new_comps, lcc_copy, m_max_possible_comps, 
                                              m_initial_volume_ch, m_depth + 1, m_sigma, m_sym_axis);
                    node_children[i] = std::shared_ptr<const Node>(new_node);
                    std::cout << "expand a new node (#comp=" << node_children[i]->num_comps() << ")\n";
                }
                else
                {
                    std::cout << "no new cut\n";
                    node_children[i] = nullptr; 
                }
            }

            std::cout << " expand done\n";
            std::cout << "num children=" << node_children.size() << '\n';
            node_children.erase(std::remove(node_children.begin(), node_children.end(), nullptr), node_children.end());
            // std::cout << "erased null children\n";
            std::cout << "num children=" << node_children.size() << '\n';
            // std::vector<int> idxs(node_children.size());
            // std::iota(idxs.begin(), idxs.end(), 0);
            // std::sort(idxs.begin(), idxs.end(),
            //     [&node_children](int a, int b) -> bool
            //     {
            //         if (node_children[a] != nullptr && node_children[b] != nullptr)
            //             return *node_children[a] < *node_children[b];
            //         if (node_children[a] == nullptr && node_children[b] != nullptr)
            //             return false;
            //         if (node_children[a] != nullptr && node_children[b] == nullptr)
            //             return true;
            //     });
            std::sort(node_children.begin(), node_children.end(),
                [](std::shared_ptr<const Node> a, std::shared_ptr<const Node> b) -> bool
                {
                    if (a != nullptr && b != nullptr)
                        return *a < *b;
                    if (a == nullptr && b != nullptr)
                        return false;
                    if (a != nullptr && b == nullptr)
                        return true;
                    if (a == nullptr && b == nullptr)
                        return false;
                });
            // std::cout << "sort " << node_children.size() << " children\n";
            if (Global_var::max_children > 0)
            {
                size_t max_children = Global_var::max_children < node_children.size() ? Global_var::max_children : node_children.size();
                node_children.erase(node_children.begin() + max_children, node_children.end());
            }
            // else Global_var::max_children < 0, keep all children
            
            // std::cout << "final num children=" << node_children.size() << '\n';
            children.insert(children.begin(), node_children.begin(), node_children.end());
        }

        Point transform(Point p)
        {
            auto &[rho, theta] = m_sym_axis;
            float a1 = 1 - 2 * std::pow(std::cos(theta), 2),
                  b1 = -2 * std::sin(theta) * std::cos(theta),
                  d1 = 2 * rho * std::cos(theta);
            float x = a1 * CGAL::to_double(p.x()) + b1 * CGAL::to_double(p.y()) + d1;
            
            float a2 = -2 * std::sin(theta) * std::cos(theta),
                b2 = 1 - 2 * std::pow(std::sin(theta), 2),
                d2 = 2 * rho * std::sin(theta);
            float y = a2 * CGAL::to_double(p.x()) + b2 * CGAL::to_double(p.y()) + d2;

            return Point(x, y, p.z());
        }

        void update_IoS()
        {
            assert (m_is_convex.size() == m_comps.size());
            m_IoS.clear();
            m_IoS.resize(m_comps.size());

            
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (!m_is_convex[i]) 
                    m_IoS[i] = 0;
                else
                {
                    Dart_const_descriptor cell_dart;
                    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
                        it=m_lcc.one_dart_per_cell<3>().begin(),
                        itend=m_lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
                    {
                        if (m_lcc.info<3>(it).first == m_comps[i]->m_in_cell)   
                        {    
                            cell_dart = it;
                            break; 
                        }
                    }
                    
                    m_IoS[i] = IoS(m_lcc, cell_dart, Global_var::sampled_unit_volume * 20, 
                            int(Global_var::min_points/20), int(Global_var::max_points/20));
                }
                    
            }
        }

        void update_degree()
        {
            m_degrees.clear();
            m_degrees.resize(m_comps.size());
            for (size_t i = 0; i < m_comps.size(); i++)
                m_degrees[i] = num_nbh_cell3(m_lcc, m_comps[i]->m_in_cell);
        }

        float get_ratio_comp_degree_2_stable() const
        {
            float ratio = 0;
            for (size_t i = 0; i < m_degrees.size(); i++)
            {
                if (m_degrees[i] == 2 && m_is_convex[i])
                {
                    auto nbhs = nbh_cell3(m_lcc, m_comps[i]->m_in_cell);
                    assert (nbhs.size() == 2);

                    int cid1 = get_compid_from_vid(*nbhs.begin()),
                        cid2 = get_compid_from_vid(*std::next(nbhs.begin()));
                    bool is_convex1, is_convex2;
                    if (cid1 >= 0) is_convex1 = m_is_convex[cid1];
                    else
                    {
                        Component* comp = new Component(*nbhs.begin(), m_lcc);
                        is_convex1 = comp->num_cc_edges() == 0;
                        
                    }   
                    if (cid2 >= 0) is_convex2 = m_is_convex[cid2];
                    else
                    {
                        Component* comp = new Component(*std::next(nbhs.begin()), m_lcc);
                        is_convex2 = comp->num_cc_edges() == 0;
                    }

                    if (is_convex1 && is_convex2)
                    {
                        ratio += m_volumes[i] / m_sum_volume;
                        #if DEBUG_CENTRAL==1
                        std::cout << "find one stable degree 2:" << 
                         m_volumes[i] << "/" << m_sum_volume 
                         << "=" << m_volumes[i] / m_sum_volume << '\n';
                        #endif
                    }
                }
            }
            #if DEBUG_CENTRAL==1
            std::cout << "volume ratio of stable degree 2 =" << ratio << '\n';
            #endif
            return ratio;
        }

        float get_ratio_comp_more_than_degree_2_stable() const
        {
            float ratio = 0;
            for (size_t i = 0; i < m_degrees.size(); i++)
            {
                if (m_degrees[i] > 2 && m_is_convex[i])
                {
                    ratio += (m_degrees[i] - 2) * m_volumes[i] / m_sum_volume;
                    #if DEBUG_DECENTRAL==1
                    std::cout << "find one stable degree >2:" << 
                        m_volumes[i] << "/" << m_sum_volume 
                        << "=" << m_volumes[i] / m_sum_volume << '\n';
                    #endif
                }
            }
            #if DEBUG_DECENTRAL==1
            std::cout << "volume ratio of stable degree >2 =" << ratio << '\n';
            #endif 
            return ratio;
        }

        int get_compid_from_vid(const Volume_idx vid) const 
        {
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (m_comps[i]->m_in_cell == vid)
                    return i;
            }
            std::cerr << "cannot find vid=" << vid 
                     << " from the node's components!\n";
            return -1;
        }

        float max_volume() const 
        {
            // if (m_volumes.size() != m_comps.size())
            // update_volume();

            return *std::max_element(m_volumes.begin(), m_volumes.end());
        }

        int max_vol_concave_comp() const
        {
            // if (m_is_convex.size() != m_comps.size())
            //     update_convexity();
            
            int max_concave = 0;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (m_is_convex[i]) continue;
                if (m_volumes[i] > m_volumes[max_concave])
                    max_concave = i;
            }
            return max_concave;
        }

        std::vector<size_t> sort_comps_to_cut() const
        {
            std::vector<size_t> comp_to_cut;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (m_is_convex[i]) continue;
                comp_to_cut.push_back(i);
            }
            std::sort(comp_to_cut.begin(), comp_to_cut.end(),
                [&](size_t a, size_t b) {
                    return m_volumes[a] > m_volumes[b];
                });
            return comp_to_cut;
        }

        int max_possible_comps() const {return m_max_possible_comps;};

        int remaining_ccv_edges() const {
            int ccv_edges = 0;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                // if (m_is_convex[i]) continue;
                // ccv_edges += m_comps[i]->num_cc_edges() + m_comps[i]->num_filtered_out_cc_edges();
                ccv_edges += m_comps[i]->num_all_cc_edges();
            }

            return ccv_edges;
        };

        float convex_sym_iou() const
        {
            float sum = 0;
            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (!m_is_convex[i]) continue;

                #if DEBUG_SYMMETRY==1
                std::cout << "IoS[" << i << "]=" << m_IoS[i] 
                          << " vol=" << m_volumes[i]/m_sum_volume 
                          << " add=" << (1 - m_IoS[i]) * ( m_volumes[i] / m_sum_volume)
                          << '\n';
                #endif

                sum += (1 - m_IoS[i]) * ( m_volumes[i] / m_sum_volume);
            }
            // sum /= m_sum_volume;
            return sum;
        }

        float sum_volume() const
        {
            // if (m_sum_volume == 0)
            //     update_volume();
            return m_sum_volume;
        };

        int max_degree() const 
        {
            return *std::max_element(m_degrees.begin(), m_degrees.end());
        };

        float accum_volume_loss() const
        {
            std::cout << "accum_volume_loss\n";

            float volume_loss = 0;
            std::vector<int> max_idx(m_comps.size());
            std::iota(max_idx.begin(), max_idx.end(), 0);
            std::sort(max_idx.begin(), max_idx.end(), 
                [&](int a, int b) {
                    return m_volumes[a] > m_volumes[b];
                });
            std::cout << "sorted\n";

            for (size_t i = 0; i < m_comps.size(); i++)
            {
                if (i > 0 && !m_is_convex[max_idx[i-1]]) break;
                volume_loss += m_sum_volume - m_volumes[max_idx[i]];
                for (size_t j = 0; j < i; j++)
                    volume_loss -= m_volumes[max_idx[j]];
            }
            return volume_loss;
        }

        int num_end_comps() const
        {
            int num = 0;
            for (size_t i = 0; i < m_comps.size(); i++)
            {   
                #if DEBUG_CENTRAL==1
                std::cout << "degree[" << i << "]=" << m_degrees[i] << '\n';
                #endif
                if (m_degrees[i] == 1)
                    num ++;
            }
            return num;
        }

        float delta_ch() const { return m_sum_delta_ch; }
        float initial_delta_ch() const {return m_initial_delta_ch;}
        float reduced_delta_ch() const {return m_initial_delta_ch - m_sum_delta_ch;}

        void save(std::string name) const
        {
            std::string off_file = name + ".off";
            std::string sta_file = name + ".txt";
            LCC_3 save_lcc = m_lcc;
            save_lcc.onmerge_functor<2>()=Merge_functor_2();

            // clean_lcc_edges(save_lcc);
            // #if VIS==1
            // CGAL::draw(save_lcc, "save-clean edges");
            // #endif
            write_lcc_as_off(save_lcc, off_file, true, true, sta_file);
            std::ofstream sta_o(sta_file, std::ios_base::app);

            sta_o << "sigma=" << m_sigma << '\n';
            sta_o << "cost=" << this->cost() << "(g=" << this->g() << " h=" << this->h() << ")\n";
            sta_o << "num comps=" << m_comps.size() << ", max possible comps=" << m_max_possible_comps << '\n';
            sta_o << "delta_ch=" << this->delta_ch() << ", reduced_delta_ch=" << this->reduced_delta_ch() << ", initial_delta_ch=" << this->initial_delta_ch() << '\n';
            sta_o << "remaining_cc_edges=" << this->remaining_ccv_edges() << '\n';
            sta_o << "max_volume=" << this->max_volume() << " sum_volume=" << this->sum_volume() << '\n';
            // TODO: write down convex sym iou even though this decomposition is not for symmetry
            sta_o << "get_ratio_comp_degree_2_stable=" << this->get_ratio_comp_degree_2_stable() << '\n';
            sta_o << "get_ratio_comp_more_than_degree_2_stable=" << this->get_ratio_comp_more_than_degree_2_stable() << '\n';
            sta_o << "max_degree=" << this->max_degree() << " num_end_comps=" << this->num_end_comps() << '\n';
            if (m_sym_axis.first >= 0)
                sta_o << "convex_sym_iou=" << this->convex_sym_iou() << '\n';
        }
    
    public:
        LCC_3 m_lcc;
        std::vector<Component*> m_comps;
        std::vector<bool> m_is_convex;
        std::vector<float> m_volumes;
        std::vector<float> m_IoS;
        std::vector<int> m_degrees;
        std::pair<float, float> m_sym_axis;
        float m_sum_volume;
        int m_max_possible_comps;
        float m_initial_volume_ch;
        float m_initial_delta_ch;
        float m_sum_delta_ch;
        size_t m_depth;
        float m_sigma = 1.0;
        
};