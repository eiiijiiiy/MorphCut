#pragma once

#include "util.hpp"
#include "Edge.hpp"
#include "global_variables.hpp"


class Component
{
public:

    static float cvx_thresh_;
    static int max_try_;
    static std::vector<Volume_idx> finalized_comps_;
    static Cut_mode m_s_cut_mode;
    static bool m_with_bisector;
    static float m_min_cut_factor;

    Component(
        Volume_idx in_cell, const LCC_3 &lcc, bool with_additional=true): 
        m_in_cell(in_cell), m_with_additional(with_additional)
    {
        m_volume = volume_of_cell3(lcc, in_cell);
        std::cout << "volume=" << m_volume << '\n';
        m_volume_of_ch = convex_hull_volume_of_cell3(lcc, in_cell, m_ch_mesh);
        std::cout << "volume of ch=" << m_volume_of_ch << '\n';
        find_concave_edges(Global_var::concave_thresh * M_PI, lcc);
        std::cout << "num_cc_edges=" << m_cc_edges.size() << '\n';
        if (m_cc_edges.m_edges.size() > 0)
            find_potential_cuts(lcc, with_additional);
        std::cout << "[New component] #cc_edge=" << m_cc_edges.size() 
                  << "(full=" << m_num_all_cc_edges << ")"
                  << ", #potential_cuts=" << m_potential_cuts.size() << '\n';

    }

    Component(const Component *other):
        m_in_cell(other->m_in_cell), 
        m_volume(other->m_volume),
        m_volume_of_ch(other->m_volume_of_ch),
        m_cc_edges(other->m_cc_edges),
        m_cc_edge_plane_areas(other->m_cc_edge_plane_areas),
        m_potential_cuts(other->m_potential_cuts),
        m_with_additional(other->m_with_additional){}

    ~Component(){};

    float volume() const { return m_volume; };
    float volume_ch() const { return m_volume_of_ch; };
    float convexity() const
    {
        if (m_volume_of_ch < 0.9 * m_volume) 
        {
            std::cerr << m_volume_of_ch << "," << m_volume << " m_volume_of_ch < m_volume, set convexity to 0\n";
            return 0;
        }
        return m_volume / m_volume_of_ch;
    };

    int num_cc_edges() const { return m_cc_edges.size(); };
    int num_filtered_out_cc_edges() const {return m_num_filtered_out_cc_edges;};
    int num_all_cc_edges() const {return m_num_all_cc_edges;}
    float delta_ch() const {
        // if (convexity() >= cvx_thresh_ || m_potential_cuts.size() == 0) return 0;
        if (convexity() < 1e-3) return m_volume;
        else return m_volume_of_ch - m_volume;
    };

    int num_potential_cuts() const { return m_potential_cuts.size(); };
    void clear_cuts() { m_potential_cuts.clear(); };

    const Plane* get_a_random_cut_and_remove_it()
    {
        // find a random cuts
        size_t num_cuts = m_potential_cuts.size();
        if (num_cuts == 0) return nullptr;
        size_t idx = rand() % num_cuts;
        const Plane* plane = m_potential_cuts[idx];
        m_potential_cuts.erase(m_potential_cuts.begin() + idx);

        return plane;
    };

    const Plane* get_first_cut_and_remove_it()
    {
        size_t num_cuts = m_potential_cuts.size();
        if (num_cuts == 0) return nullptr;
        const Plane* plane = m_potential_cuts[0];
        m_potential_cuts.erase(m_potential_cuts.begin());

        return plane;
    };

    const Plane* get_cut(int idx) const
    {
        assert (idx < m_potential_cuts.size());
        return m_potential_cuts[idx];
    };

    void remove_cut(int idx)
    {
        assert (idx < m_potential_cuts.size());
        m_potential_cuts.erase(m_potential_cuts.begin() + idx);
    };

    void remove_cut(const Plane* plane)
    {
        auto it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), plane);
        if (it != m_potential_cuts.end())
            m_potential_cuts.erase(it);
    };

    int cut(const Plane* plane, LCC_3 &lcc, std::vector<Component*> &sub_comps) const
    {
        // std::cout << "cutting component by plane\n";
        std::vector<Dart_descriptor> pos_cells, neg_cells;

        Dart_descriptor cell_dart;
        for (LCC_3::One_dart_per_cell_range<3>::iterator
            it=lcc.one_dart_per_cell<3>().begin(),
            itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
        {
            if (lcc.info<3>(it).first == m_in_cell)   
            {    
                cell_dart = it;
                break; 
            }
        }
        
        int result = cut_cell3_by_plane(
            lcc, cell_dart, plane, pos_cells, neg_cells);
        if (result != EXIT_SUCCESS)
        {
            std::cerr << "cut error\n";
            // remove_cut(plane);
            return result;
        }

        std::vector<Volume_idx> pos_vids, neg_vids;
        for (auto d : pos_cells)
            pos_vids.push_back(lcc.info<3>(d).first);
        for (auto d : neg_cells)
            neg_vids.push_back(lcc.info<3>(d).first);

        create_components(
            pos_vids, lcc, std::back_inserter(sub_comps));
        create_components(
            neg_vids, lcc, std::back_inserter(sub_comps));

        return EXIT_SUCCESS;
    };

    void find_concave_edges(const float &radian_thresh, const LCC_3 &lcc)
    {
        // std::cout << "finding concave edges ...\n";
        m_cc_edges.clear();
        // collect all edges of the components
        using Edge_idx = int;
        
        Dart_const_descriptor cell_dart;
        for (LCC_3::One_dart_per_cell_range<3>::const_iterator
            it=lcc.one_dart_per_cell<3>().begin(),
            itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
        {
            if (lcc.info<3>(it).first == m_in_cell)   
            {
                cell_dart = it;
                break; 
            }
        }

        std::map<Edge_idx, Dart_const_descriptor> edges;
        for (LCC_3::One_dart_per_incident_cell_range<1, 3>::const_iterator
            it=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).begin(),
            itend=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).end(); it!=itend; ++it)
            edges.insert({lcc.info<1>(it), it});
        
        #if DEBUG_CONCAVE_EDGE == 1
        std::cout << "num edges: " << edges.size() << '\n';
        #endif
        
        for (auto [eid, e_dart] : edges)
        {
            bool flip1, flip2;
            float inward_radian = inner_dihedral_angle(lcc, e_dart, cell_dart, flip1, flip2);
            if (inward_radian < 1e-3) continue;
            
            Dart_const_descriptor other_face_dart = lcc.beta<2>(e_dart);
            const Plane *plane1 = lcc.info<2>(e_dart).second,
                        *plane2 = lcc.info<2>(other_face_dart).second;
            if (plane1 == plane2) continue;
            std::pair<const Plane*, const Plane*> plane_pair = plane1 < plane2 ? 
                std::make_pair(plane1, plane2) : std::make_pair(plane2, plane1);
            std::pair<bool, bool> flip_pair = plane1 < plane2 ? 
                std::make_pair(flip1, flip2) : std::make_pair(flip2, flip1);

            m_cc_edges.insert(plane_pair, flip_pair, eid, e_dart, inward_radian); 
            // if (inward_radian > radian_thresh)
            // {
            //     m_cc_edge_plane_areas.insert({plane1, 0});
            //     m_cc_edge_plane_areas.insert({plane2, 0});

            //     #if DEBUG_CONCAVE_EDGE == 1
            //     bool mark1 = false, mark2 = false;
            //     size_type amark = lcc.get_new_mark();
            //     for (LCC_3::Dart_of_cell_range<2>::const_iterator
            //         it=lcc.darts_of_cell<2>(e_dart).begin(),
            //         itend=lcc.darts_of_cell<2>(e_dart).end(); it!=itend; ++it)
            //     {
            //         lcc.mark(it, amark);
            //         mark1 = true;
            //     }
            //     for (LCC_3::Dart_of_cell_range<2>::const_iterator
            //         it=lcc.darts_of_cell<2>(other_face_dart).begin(),
            //         itend=lcc.darts_of_cell<2>(other_face_dart).end(); it!=itend; ++it)
            //     {
            //         lcc.mark(it, amark);
            //         mark2 = true;
            //     }
            //     std::cout << "mark1=" << mark1 << ", mark2=" << mark2 << '\n';
            //     LCC_gs_options gso(amark);
            //     std::string edge_title = "edge-inward_radian=" + std::to_string(inward_radian);
            //     CGAL::draw(lcc, gso, edge_title.c_str());
            //     lcc.free_mark(amark);
            //     #endif
            // }
        }
        m_cc_edges.update_edge_ends(lcc);
        m_cc_edges.merge_edges(lcc, m_in_cell, Global_var::dist_thresh);
        m_cc_edges.remove_convex_edges(radian_thresh);

        Tree ch_tree = Tree(faces(m_ch_mesh).first, faces(m_ch_mesh).second, m_ch_mesh);
        m_cc_edges.update_weights(ch_tree);
        m_num_all_cc_edges = m_cc_edges.size();
        m_cc_edges.filter_edges(m_min_cut_factor);
        m_num_filtered_out_cc_edges = m_num_all_cc_edges - m_cc_edges.size();
        m_cc_edges.sort();

        #if DEBUG_CONCAVE == 1
        std::cout << "num cc edges before filtering: " << m_cc_edges.size() << '\n';
        #endif        
        
        #if DEBUG_CONCAVE_EDGE == 1
        for (auto edge : m_cc_edges.m_edges)
        {
            LCC_3 draw_lcc = lcc;
            Volume_idx cell_idx = lcc.info<3>(cell_dart).first;
            Dart_descriptor cell_dart_of_draw;
            for (LCC_3::One_dart_per_cell_range<3>::iterator
                it=draw_lcc.one_dart_per_cell<3>().begin(),
                itend=draw_lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
            {
                if (draw_lcc.info<3>(it).first == cell_idx)   
                {
                    cell_dart_of_draw = it;
                    break; 
                }
            }

            auto plane1 = edge.m_plane1, plane2 = edge.m_plane2;
            size_type amark = draw_lcc.get_new_mark();
            for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator
                it=draw_lcc.one_dart_per_incident_cell<2, 3>(cell_dart_of_draw).begin(),
                itend=draw_lcc.one_dart_per_incident_cell<2, 3>(cell_dart_of_draw).end(); it!=itend; ++it)
            {
                if (draw_lcc.info<2>(it).second == plane1 || draw_lcc.info<2>(it).second == plane2)
                {
                    for (LCC_3::Dart_of_cell_range<2>::const_iterator
                        itt=draw_lcc.darts_of_cell<2>(it).begin(),
                        ittend=draw_lcc.darts_of_cell<2>(it).end(); itt!=ittend; ++itt)
                        draw_lcc.mark(itt, amark);
                }
            }
            LCC_gs_options gso(amark);
            std::string edge_title = "edge-after-filtering" + std::to_string(edge.m_weight)
                                         + "-" + std::to_string(edge.m_cm_edges.size())
                                         + "-" + std::to_string(edge.m_distance_to_ch)
                                         + "-" + std::to_string(edge.m_length);
            CGAL::draw(draw_lcc, gso, edge_title.c_str());
            draw_lcc.free_mark(amark);
        }
        #endif
        
        #if DEBUG_CONCAVE == 1
        std::cout << "num cc edges after filtering: " << m_cc_edges.size() << '\n';
        #endif
    };    

    void find_potential_cuts(const LCC_3 &lcc, bool with_additional)
    {
        // 1. add all surface planes first
        m_potential_cuts.clear();
        for (auto edge : m_cc_edges.m_edges)
        {
            auto plane1 = edge.m_plane1, plane2 = edge.m_plane2;
            if (std::find(
                m_potential_cuts.begin(), m_potential_cuts.end(), plane1) == m_potential_cuts.end())
            {
                m_potential_cuts.push_back(plane1);
                m_potential_cut_weights[plane1] = edge.m_weight;
            }
            else
                m_potential_cut_weights[plane1] += edge.m_weight;
                
            if (std::find(
                m_potential_cuts.begin(), m_potential_cuts.end(), plane2) == m_potential_cuts.end())
            {
                m_potential_cuts.push_back(plane2);
                m_potential_cut_weights[plane2] = edge.m_weight;
            }
            else
                m_potential_cut_weights[plane2] += edge.m_weight;
        }
        std::vector<const Plane*> concave_planes;
        for (const auto &[plane, area] : m_cc_edge_plane_areas)
            concave_planes.push_back(plane);

        size_t num_surface_cuts = m_potential_cuts.size();
        // add additional plane in terms of the intersection of two planes
        size_t num_horizontal_inter = 0, num_incline_inter = 0;
        if (with_additional) {
        for (auto edge : m_cc_edges.m_edges)
        {
            auto plane1 = edge.m_plane1, plane2 = edge.m_plane2;
            float pitch1 = get_pitch(plane1), pitch2 = get_pitch(plane2);
            // float area1 = m_cc_edge_plane_areas.at(plane1), area2 = m_cc_edge_plane_areas.at(plane2);
            float a1 = CGAL::to_double(plane1->a()), b1 = CGAL::to_double(plane1->b()), 
                  c1 = CGAL::to_double(plane1->c()), d1 = CGAL::to_double(plane1->d());
            float a2 = CGAL::to_double(plane2->a()), b2 = CGAL::to_double(plane2->b()), 
                  c2 = CGAL::to_double(plane2->c()), d2 = CGAL::to_double(plane2->d());
            Eigen::Vector3f dir_inter(b1*c2 - c1*b2, c1*a2 - a1*c2, a1*b2 - b1*a2);
            dir_inter.normalize(); 

            // prepare the two bisector planes and one additional veritcal planes 
            // that most cases will use
            const Plane *bisector1 = normalize_plane(new Plane(CGAL::bisector(*plane1, *plane2)));
            const Plane* bisector2 = normalize_plane(new Plane(CGAL::bisector(plane1->opposite(), *plane2)));
            auto [bct1_exist, bct1] = is_duplicate_plane(
                Global_var::box, bisector1, m_potential_cuts, 
                Global_var::angle_thresh, Global_var::dist_thresh);
            bisector1 = bct1_exist ? bct1 : bisector1;
            auto [bct2_exist, bct2] = is_duplicate_plane(
                Global_var::box, bisector2, m_potential_cuts, 
                Global_var::angle_thresh, Global_var::dist_thresh);
            bisector2 = bct2_exist ? bct2 : bisector2;

            const Plane *vertical_plane = nullptr;
            if (!is_vertical(pitch1) && !is_vertical(pitch2))
            {
                // calculate the additional vertical plane
                float z = (b2 * d1 - d2 * b1) / (c2 * b1 - c1 * b2);
                float y = (-d1 - c1 * z) / b1;
                Eigen::Vector3f pt_on_line(0, y, z);
                Eigen::Vector3f another_pt_on_line = pt_on_line + dir_inter;
                Eigen::Vector3f vertical_pt = pt_on_line;
                vertical_pt(2) -= 100.0;

                Point pt1(pt_on_line(0), pt_on_line(1), pt_on_line(2)), 
                      pt2(another_pt_on_line(0), another_pt_on_line(1), another_pt_on_line(2)),
                      vpt(vertical_pt(0), vertical_pt(1), vertical_pt(2));
                vertical_plane = normalize_plane(new Plane(pt1, pt2, vpt));
                
                auto[vp_exist, vp] = is_duplicate_plane(
                        Global_var::box, vertical_plane, m_potential_cuts,
                        Global_var::angle_thresh, Global_var::dist_thresh);
                vertical_plane = vp_exist ? vp : vertical_plane;
            }

            if (std::abs(dir_inter(2)) > 1e-1) // non horizontal intersection
            {
                num_incline_inter++;
                if (is_vertical(pitch1) && is_vertical(pitch2))
                {
                    if (!m_with_bisector) continue;
                    bool flip1 = edge.m_flip1, flip2 = edge.m_flip2;
                    float a1 = CGAL::to_double(plane1->a()), b1 = CGAL::to_double(plane1->b()), 
                          c1 = CGAL::to_double(plane1->c()), d1 = CGAL::to_double(plane1->d());
                    float a2 = CGAL::to_double(plane2->a()), b2 = CGAL::to_double(plane2->b()),
                          c2 = CGAL::to_double(plane2->c()), d2 = CGAL::to_double(plane2->d());   
                    if (flip1) { a1 = -a1; b1 = -b1; c1 = -c1; d1 = -d1; }
                    if (flip2) { a2 = -a2; b2 = -b2; c2 = -c2; d2 = -d2; }
                    float ab1 = CGAL::to_double(bisector1->a()), bb1 = CGAL::to_double(bisector1->b()), 
                          cb1 = CGAL::to_double(bisector1->c()), db1 = CGAL::to_double(bisector1->d());
                    float ab2 = CGAL::to_double(bisector2->a()), bb2 = CGAL::to_double(bisector2->b()),
                          cb2 = CGAL::to_double(bisector2->c()), db2 = CGAL::to_double(bisector2->d());
                    Eigen::Matrix2f mat1, mat2;
                    mat1 << ab1, bb1, a1, b1;
                    mat2 << ab2, bb2, a1, b1;
                    Eigen::Vector2f vec1(-db1, -d1-2), vec2(-db2, -d1-2);
                    Eigen::Vector2f res1 = mat1.inverse() * vec1, res2 = mat2.inverse() * vec2;
                    bool same_dir_b1 = a2 * res1(0) + b2 * res1(1) + d2 < 0,
                         same_dir_b2 = a2 * res2(0) + b2 * res2(1) + d2 < 0;
                    if (same_dir_b1 == same_dir_b2)
                    {
                        std::cout << "same_dir_b1 == same_dir_b2 what happen??\n";
                        continue;
                    }
                    const Plane* picked_bisector = same_dir_b1 ? bisector1 : bisector2;

                    // get two bisectors
                    std::vector<const Plane*>::iterator 
                        it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), picked_bisector);

                    if (it == m_potential_cuts.end())
                    {
                        m_potential_cuts.push_back(picked_bisector);
                        m_potential_cut_weights[picked_bisector] = edge.m_weight;
                        
                    }
                    else // it1 != m_potential_cuts.end()
                        m_potential_cut_weights[picked_bisector] += edge.m_weight;
                }
                // else if (is_vertical(pitch1) || is_vertical(pitch2)) // no thing to do 
                else if (is_incline(pitch1) && is_incline(pitch2))
                {
                   // get one vertical plane
                    std::vector<const Plane*>::iterator 
                        it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), vertical_plane);
                    if (it == m_potential_cuts.end())
                    {
                        m_potential_cuts.push_back(vertical_plane);
                        m_potential_cut_weights[vertical_plane] = edge.m_weight;
                    }
                    else // it != m_potential_cuts.end()
                        m_potential_cut_weights[vertical_plane] += edge.m_weight;
                }
                // else : one incline and one vertical
            }
            else
            {
                num_horizontal_inter++;
                if ((is_horizontal(pitch1) && is_vertical(pitch2))
                 || (is_horizontal(pitch2) && is_vertical(pitch1)))
                {
                    // get one bisector 
                    // (TODO: we only need one bisector for a such edge, but it's hard to distinguish globally)
                    // So add both 
                    // std::cout << "add two bisectors!\n";
                    // std::vector<const Plane*>::iterator 
                    //     it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), bisector1);
                    // if (it == m_potential_cuts.end())
                    // {
                    //     m_potential_cuts.push_back(bisector1);
                    //     m_cc_edge_plane_areas.insert({bisector1, (area1 + area2) / 2});
                    // }
                    // else // it1 != m_potential_cuts.end()
                    // {
                    //     int idx = std::distance(m_potential_cuts.begin(), it);
                    //     if (idx > num_surface_cuts)
                    //         m_cc_edge_plane_areas[bisector1] += (area1 + area2) / 2;
                    // }

                    // it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), bisector2);
                    // if (it == m_potential_cuts.end())
                    // {
                    //     m_potential_cuts.push_back(bisector2);
                    //     m_cc_edge_plane_areas.insert({bisector2, (area1 + area2) / 2});
                    // }
                    // else // it2 != m_potential_cuts.end()
                    // {
                    //     int idx = std::distance(m_potential_cuts.begin(), it);
                    //     if (idx > num_surface_cuts)
                    //         m_cc_edge_plane_areas[bisector2] += (area1 + area2) / 2;
                    // }
                }
                else if (is_horizontal(pitch1) || is_horizontal(pitch2))
                {
                    std::cout << "at least one horizontal?\n";
                    // get one vertical plane
                    std::vector<const Plane*>::iterator 
                        it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), vertical_plane);
                    if (it == m_potential_cuts.end())
                    {
                        std::cout << "add one vertical plane for at least one horizontal\n";
                        m_potential_cuts.push_back(vertical_plane);
                        m_potential_cut_weights[vertical_plane] = edge.m_weight;
                    }
                    else // it != m_potential_cuts.end()
                        m_potential_cut_weights[vertical_plane] += edge.m_weight;
                    
                }
                else if (is_incline(pitch1) && is_incline(pitch2))// two incline planes intersecting at a horizontal line
                {
                    // std::cout << "horizontal intersection of two incline planes\n";
                    float z = (d2 * b1 - d1 * b2) / (c1 * b2 - c2 * b1);
                    const Plane* horizontal_plane = new Plane(0, 0, 1, -z); // 0 * x + 0 * y + 1 * z - z = 0
                    auto [hp_exist, hp] = is_duplicate_plane(
                        Global_var::box, horizontal_plane, m_potential_cuts, 
                        Global_var::angle_thresh, Global_var::dist_thresh);
                    horizontal_plane = hp_exist ? hp : horizontal_plane;
                    std::vector<const Plane*>::iterator 
                        it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), horizontal_plane);
                    if (it == m_potential_cuts.end())
                    {
                        m_potential_cuts.push_back(horizontal_plane);
                        m_potential_cut_weights[horizontal_plane] = edge.m_weight;
                    }
                    else // it != m_potential_cuts.end()
                        m_potential_cut_weights[horizontal_plane] += edge.m_weight;

                    // get one vertical plane
                    it = std::find(m_potential_cuts.begin(), m_potential_cuts.end(), vertical_plane);
                    if (it == m_potential_cuts.end())
                    {
                        m_potential_cuts.push_back(vertical_plane);
                        m_potential_cut_weights[vertical_plane] = edge.m_weight;
                    }
                    else // it != m_potential_cuts.end()
                        m_potential_cut_weights[vertical_plane] += edge.m_weight;
                }
                // else : one incline and one vertical
            }
        }
        std::cout << "num_incline_inter=" << num_incline_inter 
                << ", num_horizontal_inter=" << num_horizontal_inter << '\n';
        }
        // sort the potential cuts by the area of the cc edge planes
        std::sort(m_potential_cuts.begin(), m_potential_cuts.end(), 
            [this](const Plane* p1, const Plane* p2)
            {
                return m_potential_cut_weights.at(p1) > m_potential_cut_weights.at(p2);
            });
        
        // filter by cut mode
        // std::cout << "m_s_cut_mode=" << m_s_cut_mode << "\n";
        std::cout << "#potential_cuts=" << m_potential_cuts.size();
        switch (m_s_cut_mode)
        {
            case Cut_mode::HV:
                m_potential_cuts.erase(
                    std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
                        [this](const Plane* p) -> bool
                        {
                            return is_incline(get_pitch(p));
                        }), m_potential_cuts.end());
                break;
            case Cut_mode::H:
                m_potential_cuts.erase(
                    std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
                        [this](const Plane* p) -> bool
                        {
                            return !is_horizontal(get_pitch(p));
                        }), m_potential_cuts.end());
                break;

            case Cut_mode::V:
                m_potential_cuts.erase(
                    std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
                        [this](const Plane* p) -> bool
                        {
                            return !is_vertical(get_pitch(p));
                        }), m_potential_cuts.end());
                break;

            case Cut_mode::HI:
                m_potential_cuts.erase(
                    std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
                        [this](const Plane* p) -> bool
                        {
                            return is_vertical(get_pitch(p));
                        }), m_potential_cuts.end());
                break;

            case Cut_mode::VI:
                m_potential_cuts.erase(
                    std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
                        [this](const Plane* p) -> bool
                        {
                            return is_horizontal(get_pitch(p));
                        }), m_potential_cuts.end());
                break;

            default:
                break;
        }
        std::cout << ", after mode selection=" << m_potential_cuts.size() << "\n";

        // suppress parallel cuts
        // size_t num_to_suppress = 0;
        // std::map<const Plane*, bool> to_suppress;
        // for (const auto &plane : m_potential_cuts)
        //     to_suppress.insert({plane, false});
        // for (size_t i = 0; i < m_potential_cuts.size(); ++i)
        // {
        //     if (to_suppress[m_potential_cuts[i]]) continue;
        //     for (size_t j = i + 1; j < m_potential_cuts.size(); ++j)
        //     {
        //         if (to_suppress[m_potential_cuts[j]]) continue;
        //         if (are_plane_parallel(m_potential_cuts[i], m_potential_cuts[j]))
        //         {
        //             to_suppress[m_potential_cuts[j]] = true;
        //             num_to_suppress++;
        //         }
        //     }
        // }
        // std::cout << "num_to_suppress(parallel)=" << num_to_suppress << '\n';

        // m_potential_cuts.erase(
        //     std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
        //         [&to_suppress](const Plane* p) -> bool
        //         {
        //             return to_suppress[p];
        //         }), m_potential_cuts.end());
        
        // num_to_suppress = 0;
        // to_suppress.clear();
        // std::map<const Plane*, std::vector<Point>> plane_intersections;
        // get_cut_intersections(lcc, plane_intersections);
        
        // for (size_t i = 0; i < m_potential_cuts.size(); ++i)
        // {
        //     if (to_suppress[m_potential_cuts[i]]) continue;
        //     for (size_t j = i + 1; j < m_potential_cuts.size(); ++j)
        //     {
        //         if (to_suppress[m_potential_cuts[j]]) continue;
        //         if (!are_cuts_intersecting_in_comp(
        //             m_potential_cuts[i], plane_intersections.at(m_potential_cuts[i]),
        //             m_potential_cuts[j], plane_intersections.at(m_potential_cuts[j])))
        //         {
        //             to_suppress[m_potential_cuts[j]] = true;
        //             num_to_suppress++;
        //         }
        //     }
        // }
        // std::cout << "num_to_suppress(non intersecting)=" << num_to_suppress << '\n';

        // m_potential_cuts.erase(
        //     std::remove_if(m_potential_cuts.begin(), m_potential_cuts.end(), 
        //         [&to_suppress](const Plane* p) -> bool
        //         {
        //             return to_suppress[p];
        //         }), m_potential_cuts.end());
    }

    void get_cut_intersections(
        const LCC_3 &lcc,
        std::map<const Plane*, std::vector<Point>> &cut_inters)
    {
        Dart_const_descriptor cell_dart;
        for (LCC_3::One_dart_per_cell_range<3>::const_iterator
            it=lcc.one_dart_per_cell<3>().begin(),
            itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
        {
            if (lcc.info<3>(it).first == m_in_cell)   
            {
                cell_dart = it;
                break; 
            }
        }

        std::unordered_map<int, Dart_const_descriptor> darts_of_cell_edges;
        for (LCC_3::One_dart_per_incident_cell_range<1, 3>::const_iterator
            it=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).begin(),
            itend=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).end(); it!=itend; ++it)
            darts_of_cell_edges[lcc.info<1>(it)] = it;
        
        for (const auto plane : m_potential_cuts)
        {
            std::vector<Point> pts;
            for (const auto [eid, e] : darts_of_cell_edges)
            {
                Dart_const_descriptor vtx0 = e;
                Dart_const_descriptor vtx1 = lcc.beta<2>(e);
                LCC_3::Point p0 = lcc.point(vtx0);
                LCC_3::Point p1 = lcc.point(vtx1);
                Segment seg(Point(p0.x(), p0.y(), p0.z()), Point(p1.x(), p1.y(), p1.z()));
                
                // 1.3. Get the intersection point
                if (CGAL::do_intersect(seg, *plane))
                {
                    CGAL::Object obj = CGAL::intersection(seg, *plane);
                    const Point* p = CGAL::object_cast<Point>(&obj);
                    const Segment* s = CGAL::object_cast<Segment>(&obj);
                    if (p)
                        pts.push_back(*p);
                    if (s)
                    {
                        pts.push_back(s->source());
                        pts.push_back(s->target());
                    }
                }
            }  
            cut_inters[plane] = pts;
        }
    }

    bool are_cuts_intersecting_in_comp(
        const Plane* plane1, const std::vector<Point> &pts1,
        const Plane* plane2, const std::vector<Point> &pts2)
    {
        // project the intersection points onto plane1
        std::vector<Point2> pts1_2d, pts2_2d;
        for (const auto pt : pts1)
            pts1_2d.push_back(plane1->to_2d(pt));
        for (const auto pt : pts2)
            pts2_2d.push_back(plane1->to_2d(pt));

        // coarse intersection: check their bounding box :
        CGAL::Bbox_2 bbox1 = CGAL::bbox_2(pts1_2d.begin(), pts1_2d.end());
        CGAL::Bbox_2 bbox2 = CGAL::bbox_2(pts2_2d.begin(), pts2_2d.end());

        if (bbox1.xmin() >= bbox2.xmax() 
         || bbox2.xmin() >= bbox1.xmax()
         || bbox1.ymin() >= bbox2.ymax() 
         || bbox2.ymin() >= bbox1.ymax())
            return false;
        else
            return true;
        // TODO: make this more precise?
    }
    
    void create_components(
        const std::vector<Volume_idx> &comps,
        LCC_3 &lcc,
        std::back_insert_iterator<std::vector<Component*>> insertIt) const
    {
        // create components for the pos components
        for (size_t i = 0; i < comps.size(); ++i)
        {   
            std::cout << "creating new components after cutting\n";
            Component* comp = new Component(comps[i], lcc);
            *(insertIt++) = comp;
        }
    };

    bool is_vertical(float pitch) {
        assert(Global_var::radian_thresh > 0 && Global_var::radian_thresh < 0.5 * M_PI);
        assert(pitch > -1e-6 && pitch < 0.5 * M_PI + 1e-6);

        if (pitch < Global_var::radian_thresh)
            return true;
        else
            return false;
    }

    bool is_horizontal(float pitch) {
        assert(Global_var::radian_thresh > 0 && Global_var::radian_thresh < 0.5 * M_PI);
        assert(pitch > -1e-6 && pitch < 0.5 * M_PI + 1e-6);
        
        if (pitch > 0.5 * M_PI - Global_var::radian_thresh)
            return true;
        else
            return false;
    }

    bool is_incline(float pitch) {
        if (!is_vertical(pitch) && !is_horizontal(pitch))
            return true;
        else
            return false;
    }

    Volume_idx m_in_cell;
    float m_volume, m_volume_of_ch;
    size_t m_num_all_cc_edges;
    size_t m_num_filtered_out_cc_edges;
    EdgeSet m_cc_edges;
    std::map<const Plane*, float> m_cc_edge_plane_areas;
    std::map<const Plane*, float> m_potential_cut_weights;
    std::vector<const Plane*> m_potential_cuts;
    bool m_with_additional;
    Surface_mesh m_ch_mesh;
};