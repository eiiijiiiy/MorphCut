#pragma once

#include <queue>
#include <set>
#include <vector>
#include "util.hpp"
#include "Node.hpp"
#include "global_variables.hpp"

// TODO: constrain the max size of priority queue

void expand(
    std::shared_ptr<const Node> root, 
    std::vector<std::shared_ptr<const Node>> &children)
{
    root->expand(children);
}

std::function<float(const Node*)> g_fewer_comp = 
    [](const Node* node) -> float 
{
    return float(node->num_comps()) / node->max_possible_comps()
           + node->reduced_delta_ch() / node->initial_delta_ch();
};

std::function<float(const Node*)> g_max_volume = 
    [](const Node* node) -> float 
{
    // std::cout << " num=" << float(node->num_comps()) / node->max_possible_comps() 
    //           << " vol=" << 1 - node->max_volume() / node->sum_volume() << '\n';
    return Global_var::mv_weight * (1 - node->max_volume() / node->sum_volume());
};

std::function<float(const Node*)> g_symmetry = 
    [] (const Node* node) -> float
{
    return Global_var::sym_weight * node->convex_sym_iou() ;
};

std::function<float(const Node*)> g_central = 
    [] (const Node* node) -> float
{
    // return node->num_comps() + weight * (node->num_comps() - node->num_end_comps());
    return Global_var::cen_weight * node->get_ratio_comp_degree_2_stable();
};

std::function<float(const Node*)> g_decentral = 
    [] (const Node* node) -> float
{
    // return node->num_comps() +  weight * (std::max(2, node->max_degree()) - 2);
    return Global_var::de_weight * node->get_ratio_comp_more_than_degree_2_stable();
};

std::function<float(const Node*)> f_general_n = 
    [] (const Node* node) -> float
{
    return float(node->remaining_ccv_edges()) / (4 * node->max_possible_comps())
           + node->delta_ch() / node->initial_delta_ch();
};

struct NodeCmp
{
    bool operator()(std::shared_ptr<const Node> lhs, std::shared_ptr<const Node> rhs) const 
    {
        return *lhs < *rhs;
    }
};

int bnb(
    const Node* root,
    std::vector<std::shared_ptr<const Node>> &result,
    std::string log_search_path)
{
    std::set<std::shared_ptr<const Node>, NodeCmp> open;
    std::set<std::shared_ptr<const Node>, NodeCmp> droped;
    std::set<std::shared_ptr<const Node>, NodeCmp> closed; // how to control the cost for LOD????
    closed.clear();
    float closed_cost = std::numeric_limits<float>::max();
    
    #if LOG_SEARCH == 1
    std::ofstream log_search(log_search_path);
    log_search << "log bnb search as follows:\n";
    // log_search.close();
    #endif

    if (root->h() == 0)
    {   
        result.push_back(std::shared_ptr<const Node>(root));
        return EXIT_SUCCESS;
    }
    else
    {
        open.insert(std::shared_ptr<const Node>(root));
        while (!open.empty())
        {
            auto top = *open.begin();
            float top_cost = top->cost();
            std::cout << "[SEARCH] open: #=" << open.size() 
                      << "(top: cost=" << top_cost << " h=" << top->h() << " g=" << top->g() << ")"
                      << ", closed: #=" << closed.size() << "(cost=" << closed_cost << ")\n";
            #if LOG_SEARCH == 1
            // std::ofstream log_search(log_search_path, std::ios::app);
            log_search << "[SEARCH] open: #=" << open.size() 
                      << "(top: cost=" << top_cost << " h=" << top->h() << " g=" << top->g() << ")"
                      << ", closed: #=" << closed.size() << "(cost=" << closed_cost << ")\n";
            #endif

            #if VIS_NODE == 1
            top->draw("top");
            #endif
            droped.insert(top);
            std::cout << "[SEARCH] droped: #=" << droped.size() << '\n';
            open.erase(open.begin());
            
            std::vector<std::shared_ptr<const Node>> children;
            expand(top, children);

            #if DEBUG_SEARCH == 1
            std::cout << "[SEARCH] expand: " << children.size() << '\n';
            #endif

            if (children.empty())
            {
                #if DEBUG_SEARCH == 1
                std::cout << "[SEARCH] - Dead top! \n";
                #endif
                if (top->cost() < closed_cost - 1e-3)
                {
                    closed.clear();
                    closed.insert(top);
                    closed_cost = (*closed.begin())->cost();
                    #if VIS_NODE == 1 
                    top->draw("dead_open");
                    #endif  

                    if (!open.empty())
                    {
                        std::set<std::shared_ptr<const Node>, NodeCmp>::iterator it = open.begin();
                        float other_cost;
                        do{
                            float other_cost = (*it)->cost();
                            #if DEBUG_SEARCH == 1
                            std::cout << "[SEARCH] node #" << std::distance(open.begin(), it) << "(cost=" << other_cost << ")\n";
                            #endif
                            if (other_cost >= closed_cost) break;
                        }while((++it) != open.end());

                        if (it != open.end())
                        {
                            #if DEBUG_SEARCH == 1
                            std::cout << "[SEARCH] erase " << std::distance(it, open.end()) << "from " << open.size() << " in open\n";
                            #endif
                            open.erase(it, open.end());
                        }
                    }
                }
                else if (top->cost() >= closed_cost - 1e-3 && top->cost() < closed_cost + 1e-3)
                    closed.insert(top);
                else
                {
                    #if DEBUG_SEARCH == 1
                    std::cout << "[SEARCH] - Cost larger than closed. Drop. \n";
                    #endif
                }
                // float top_cost = top->cost();
                // if (closed.size() > 0)
                // {
                //     // if (closed[0]->h() >0) 
                //     // {
                //     if (top_cost < closed_cost)
                //     {
                //         closed.clear();
                //         closed_cost = top_cost;
                //         #if DEBUG_SEARCH == 1
                //         std::cout << "[SEARCH] Update closed_cost=" << closed_cost << '\n';
                //         #endif
                //         if (!open.empty())
                //         {
                //             std::set<std::shared_ptr<const Node>, NodeCmp>::iterator it = open.begin();
                //             float other_cost;
                //             do{
                //                 float other_cost = (*it)->cost();
                //                 #if DEBUG_SEARCH == 1
                //                 std::cout << "[SEARCH] node #" << std::distance(open.begin(), it) << "(cost=" << other_cost << ")\n";
                //                 #endif
                //                 if (other_cost >= closed_cost) break;
                //             }while((++it) != open.end());

                //             if (it != open.end())
                //             {
                //                 #if DEBUG_SEARCH == 1
                //                 std::cout << "[SEARCH] erase " << std::distance(it, open.end()) << "from " << open.size() << " in open\n";
                //                 #endif
                //                 open.erase(it, open.end());
                //             }
                //         }
                //         closed.push_back(top);
                //     }
                //     // }
                // }
                // else
                // {
                //     closed_cost = top_cost;
                //     #if DEBUG_SEARCH == 1
                //     std::cout << "[SEARCH] Update closed_cost=" << closed_cost << '\n';
                //     #endif
                //     closed.push_back(top);
                // }
            }
            else
            {
                // delete top;
                for (auto child : children)
                {
                    float cost = child->cost();

                    #if DEBUG_SEARCH == 1
                    std::cout << "[SEARCH] cost of child: " << cost 
                            << "(g=" << child->g() << ", h=" << child->h() << ", ccv_edges=" << child->remaining_ccv_edges() <<")\n";
                    #endif

                    #if VIS_CHILD == 1
                    child->draw("child");
                    #endif
                    
                    if (droped.find(child) != droped.end())
                    {
                        #if DEBUG_SEARCH == 1
                        std::cout << "[SEARCH] - Node found in droped!. \n";
                        #endif
                        continue;
                    }

                    if (child->h() > 0)
                    {
                        if (cost < closed_cost)
                            open.insert(child);
                        else
                        {
                            #if DEBUG_SEARCH == 1
                            std::cout << "[SEARCH] - Cost larger than closed. Drop. \n";
                            #endif
                        }
                    }
                    else
                    {
                        if (cost < closed_cost - 1e-3)
                        {
                            closed_cost = cost;
                            #if DEBUG_SEARCH == 1
                            std::cout << "[SEARCH] Update closed_cost=" << closed_cost << '\n';
                            #endif
                            #if VIS_NODE == 1
                            child->draw();
                            #endif
                            // clean the closed and open priority queue
                            closed.clear();

                            if (!open.empty())
                            {
                                std::set<std::shared_ptr<const Node>, NodeCmp>::iterator it = open.begin();
                                float other_cost;
                                do{
                                    float other_cost = (*it)->cost();
                                    #if DEBUG_SEARCH == 1
                                    std::cout << "[SEARCH] node #" << std::distance(open.begin(), it) << "(cost=" << other_cost << ")\n";
                                    #endif
                                    if (other_cost >= closed_cost) break;
                                }while((++it) != open.end());

                                if (it != open.end())
                                {
                                    #if DEBUG_SEARCH == 1
                                    std::cout << "[SEARCH] erase " << std::distance(it, open.end()) << "from " << open.size() << " in open\n";
                                    #endif
                                    open.erase(it, open.end());
                                }
                            }
                            closed.insert(child);
                            // break;
                        }
                        else if (cost >= closed_cost - 1e-3 and cost < closed_cost + 1e-3)
                            closed.insert(child);
                        else
                        {
                            // delete child;
                            #if DEBUG_SEARCH == 1
                            std::cout << "[SEARCH] - Cost larger than existing closed. Drop. \n";
                            #endif
                        }
                    }
                }
            }
            
            // if (closed.size() > 0) break;

            // if(open.size() >= max_size)
            // {
            //     Node* top = *open.begin();
            //     std::vector<Node*> node_to_remove_earlier;
            //     for (auto it = --open.end(); it != std::next(open.begin(), half_size); it--)
            //     {
            //         Node* node = *it;
            //         if (node->g() > top->g() && node->h() > top->h())
            //             node_to_remove_earlier.push_back(node);
            //     }
            //     #if DEBUG_SEARCH == 1
            //     std::cout << "[Search] remove " << node_to_remove_earlier.size() << " nodes from open earlier\n";
            //     #endif
            //     for (auto node : node_to_remove_earlier)
            //         open.erase(node);
            // }
        }
        #if LOG_SEARCH == 1
        log_search << "[Search] Done!\n\t closed: " << closed.size() 
            << "(cost: " << closed_cost << ", g=" << (*closed.begin())->g() << ", h=" << (*closed.begin())->h() << ")\n";
        log_search.close();
        #endif
        std::cout << "[Search] Done!\n\t closed: " << closed.size() 
            << "(cost: " << closed_cost << ", h=" << (*closed.begin())->h() << ", g=" << (*closed.begin())->g() << ")\n";
        result.clear();
        result = std::vector<std::shared_ptr<const Node>>(closed.begin(), closed.end());
        // result.push_back(*closed.begin());
        // if (dead_open.size() > 0)
        //         std::cout << "(top: cost=" << (*dead_open.begin())->cost() 
        //                 << " h=" << (*dead_open.begin())->h() 
        //                 << " g=" << (*dead_open.begin())->g();
        // std::cout<< ")\n";

        // if (closed.size() == 0 && dead_open.size() > 0)
        // {
        //     std::cout << "[Search] Use top dead open as results\n";
        //     closed_cost = (*dead_open.begin())->cost();
        //     while (!dead_open.empty() && (*dead_open.begin())->cost() == closed_cost)
        //     {
        //         closed.push_back(*dead_open.begin());
        //         dead_open.erase(dead_open.begin());
        //     }
        // }
        return EXIT_SUCCESS;
    }
}


int bnb_lod(
    const LCC_3 &lcc,
    const std::vector<float> &min_cut_factors,
    std::vector<std::shared_ptr<const Node>> &lod_results,
    const float sigma, 
    std::string filename,
    const std::pair<float, float> &sym_axis = {-1, 0})
{
    lod_results.clear();
    int initial_cc_edges = 0;
    float initial_volume_ch = 0;
    float initial_volume = 0;
    size_t lod = 0;
    LCC_3 cur_lcc = lcc;
    std::vector<Volume_idx> cell3; 

    do{
        // float ccv_thresh = lod_threshs[lod];
        // assert (ccv_thresh > 0 && ccv_thresh < 1.0);
        // if (lod > 0) assert(lod_threshs[lod] > lod_threshs[lod-1]);
        if (lod > 0) assert(min_cut_factors[lod] < min_cut_factors[lod-1]);
        Node::cvx_thresh_ = 0.999;
        Component::cvx_thresh_ = 0.999;
        Component::m_min_cut_factor = min_cut_factors[lod];

        // collect all 3-cells
        cell3.clear();
        for (LCC_3::One_dart_per_cell_range<3>::iterator 
            it=cur_lcc.one_dart_per_cell<3>().begin(), itend=cur_lcc.one_dart_per_cell<3>().end();
            it != itend; ++it)
            cell3.push_back(cur_lcc.info<3>(it).first);

        // decompose one cell3
        bool with_additional = true;
        for (size_t i = 0; i < cell3.size(); i++)
        {
            cur_lcc.onsplit_functor<1>() = Split_functor_1(cur_lcc);    
            cur_lcc.onsplit_functor<2>() = Split_functor_2(cur_lcc);
            cur_lcc.onsplit_functor<3>() = Split_functor_3(cur_lcc);
            cur_lcc.onmerge_functor<2>() = Merge_functor_2();
            cur_lcc.onmerge_functor<3>() = Merge_functor_3();

            Volume_idx in_cell = cell3[i];
            Component* comp = new Component(in_cell, cur_lcc, with_additional);
            if (lod == 0) 
            {
                assert (cell3.size() == 1);
                initial_cc_edges = comp->num_all_cc_edges();
                initial_volume_ch = comp->volume_ch();
                initial_volume = comp->volume();
            }

            if (comp->num_potential_cuts() > 0)
            {
                size_t num_cc_edges = comp->num_all_cc_edges();
                float volume_ch = comp->volume_ch();
                const Node* comp_node = new Node({comp}, cur_lcc, num_cc_edges + 1, 
                                                volume_ch, 0, sigma, sym_axis);
                std::vector<std::shared_ptr<const Node>> closed;
                std::string log_search_path = filename + "log_search_lod_" + std::to_string(lod) + "_cell_" + std::to_string(in_cell) + ".txt";
                if (bnb(comp_node, closed, log_search_path) == EXIT_SUCCESS && closed.size() > 0)
                    cur_lcc = closed[0]->m_lcc; 
            }
        }
        merge_small_cell(cur_lcc, initial_volume, 0.1 * min_cut_factors[lod]);
        std::cout << "merged!\n";
        #if VIS_LOD == 1
        std::string name = "lod_" + std::to_string(lod);
        CGAL::draw(cur_lcc, name.c_str());
        #endif

        std::vector<Component*> comps;
        for (LCC_3::One_dart_per_cell_range<3>::iterator 
            it=cur_lcc.one_dart_per_cell<3>().begin(), itend=cur_lcc.one_dart_per_cell<3>().end();
            it != itend; ++it)
        {
            Volume_idx in_cell = cur_lcc.info<3>(it).first;
            Component* comp = new Component(in_cell, cur_lcc, false);
            comps.push_back(comp);
        }  
        const Node* result = new Node(comps, cur_lcc, initial_cc_edges+1, 
                                     initial_volume_ch, 0, sigma, sym_axis);

        std::string name = filename + "_lod_" + std::to_string(min_cut_factors[lod]);
        result->save(name);

        lod_results.push_back(std::shared_ptr<const Node>(result));
        lod++;
    }while(lod < min_cut_factors.size());
    
    std::cout << "EXIT with lod=" << lod-1 << "(min_cut_factor=" << min_cut_factors[lod-1] << ")\n";
    return EXIT_SUCCESS;
}