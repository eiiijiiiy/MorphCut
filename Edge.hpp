#pragma once

#include "util.hpp"
#include "global_variables.hpp"

class Edge{

public: 
    Edge(const Plane* plane1, const Plane* plane2, 
         const bool flip1, const bool flip2,
         const float inward_radian)
         :m_inward_radian(inward_radian)
    {
        assert(plane1 != plane2);
        m_plane1 = plane1 < plane2 ? plane1 : plane2;
        m_plane2 = plane1 < plane2 ? plane2 : plane1;
        m_flip1 = plane1 < plane2 ? flip1 : flip2;
        m_flip2 = plane1 < plane2 ? flip2 : flip1;
        update_parameter_of_intersection();
        update_proj_from_origin();
    }

    Edge(
        const std::pair<const Plane*, const Plane*> &edge,
        const std::pair<bool, bool> &flip,
        const float inward_radian,
        const Eigen::Vector3f &inter_direction,
        const Eigen::Vector3f &inter_point,
        const float m_proj_from_origin,
        const float m_distance_from_origin,
        const Point &end1, const Point &end2)
        :m_plane1(edge.first), m_plane2(edge.second),
         m_flip1(flip.first), m_flip2(flip.second),
         m_inward_radian(inward_radian),
         m_inter_direction(inter_direction),
         m_inter_point(inter_point),
         m_proj_from_origin(m_proj_from_origin),
         m_distance_from_origin(m_distance_from_origin)
        {
            update_two_ends(end1, end2);
        }
    
    Edge(
        const std::pair<const Plane*, const Plane*> &edge,
        const std::pair<bool, bool> &flip,
        const float inward_radian,
        const Eigen::Vector3f &inter_direction,
        const Eigen::Vector3f &inter_point,
        const float m_proj_from_origin,
        const float m_distance_from_origin,
        float end_min, float end_max)
        :m_plane1(edge.first), m_plane2(edge.second),
         m_flip1(flip.first), m_flip2(flip.second), 
        m_inward_radian(inward_radian),
         m_inter_direction(inter_direction),
         m_inter_point(inter_point),
         m_proj_from_origin(m_proj_from_origin),
         m_distance_from_origin(m_distance_from_origin)
         {
            update_two_ends(end_min, end_max);
         }

    Edge (
        const std::pair<const Plane*, const Plane*> &edge, 
        const std::pair<bool, bool> &flip,
        const float inward_radian)
        :m_inward_radian(inward_radian)
    {
        m_plane1 = edge.first;
        m_plane2 = edge.second;
        m_flip1 = flip.first;
        m_flip2 = flip.second;
        update_parameter_of_intersection();
        update_proj_from_origin();
    }

    void insert_cm_edge(const Edge_idx edge, Dart_const_descriptor edge_dart)
    {
        m_cm_edges[edge] = edge_dart;
    }

    bool operator < (const Edge &other) const
    {
        if (m_plane1 != other.m_plane1) return m_plane1 < other.m_plane1;
        else return m_plane2 < other.m_plane2;
    }

    int merge(
        const LCC_3& lcc, Volume_idx cell_idx,
        const Edge &other, const float dist_thresh, std::vector<Edge> &new_edges) const
    {
        if (std::abs(m_inter_direction.dot(other.m_inter_direction)) < 0.95) return 0;
        new_edges.clear();
        if ((m_plane1 == other.m_plane1 && m_flip1 == other.m_flip1) 
            || (m_plane2 == other.m_plane1 && m_flip2 == other.m_flip1)
            || (m_plane1 == other.m_plane2 && m_flip1 == other.m_flip2)
            || (m_plane2 == other.m_plane2 && m_flip2 == other.m_flip2))
        {
            if (m_plane1 == other.m_plane1 && m_plane2 == other.m_plane2) return 0;

            float dist = ((m_inter_point + m_proj_from_origin * m_inter_direction) - 
                        (other.m_inter_point + other.m_proj_from_origin * other.m_inter_direction)).norm();
            // std::cout <<  "dist of two parallel intersections =" << dist << "\n";
            if (dist > dist_thresh) return 0;

            const Plane* new_plane1 = (m_plane1 == other.m_plane1 || m_plane1 == other.m_plane2) ? m_plane2 : m_plane1;
            bool new_flip1 = (m_plane1 == other.m_plane1 || m_plane1 == other.m_plane2) ? m_flip2 : m_flip1;
            const Plane* new_plane2 = (other.m_plane1 == m_plane1 || other.m_plane1 == m_plane2) ? other.m_plane2 : other.m_plane1;
            bool new_flip2 = (other.m_plane1 == m_plane1 || other.m_plane1 == m_plane2) ? other.m_flip2 : other.m_flip1;
            float new_inward_radian = (m_inward_radian + other.m_inward_radian) - M_PI;

            if (new_plane1 == new_plane2) return 2;

            // #if DEBUG_CONCAVE_EDGE == 1
            // LCC_3 draw_lcc = lcc;
            // Dart_descriptor cell_dart_of_draw;
            // for (LCC_3::One_dart_per_cell_range<3>::iterator
            //     it=draw_lcc.one_dart_per_cell<3>().begin(),
            //     itend=draw_lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
            // {
            //     if (draw_lcc.info<3>(it).first == cell_idx)   
            //     {
            //         cell_dart_of_draw = it;
            //         break; 
            //     }
            // }
            
            // size_type amark = draw_lcc.get_new_mark();
            // for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator
            //     it=draw_lcc.one_dart_per_incident_cell<2, 3>(cell_dart_of_draw).begin(),
            //     itend=draw_lcc.one_dart_per_incident_cell<2, 3>(cell_dart_of_draw).end(); it!=itend; ++it)
            // {
            //     if (draw_lcc.info<2>(it).second == m_plane1 
            //     || draw_lcc.info<2>(it).second == m_plane2
            //     || draw_lcc.info<2>(it).second == other.m_plane1
            //     || draw_lcc.info<2>(it).second == other.m_plane2)
            //     {
            //         for (LCC_3::Dart_of_cell_range<2>::const_iterator
            //             itt=draw_lcc.darts_of_cell<2>(it).begin(),
            //             ittend=draw_lcc.darts_of_cell<2>(it).end(); itt!=ittend; ++itt)
            //             draw_lcc.mark(itt, amark);
            //     }
            // }
            // LCC_gs_options gso(amark);
            // float dot = m_inter_direction.dot(other.m_inter_direction);
            // std::string edge_title = "coplanar edges-" 
            //                           + std::to_string(dot) +"-" + std::to_string(dist) + "-" + std::to_string(new_inward_radian);
            // CGAL::draw(draw_lcc, gso, edge_title.c_str());
            // draw_lcc.free_mark(amark);
            // #endif

            Edge merged_edge =  Edge(new_plane1, new_plane2, new_flip1, new_flip2, new_inward_radian);
            if (merged_edge.m_inter_direction.norm() == 0) return 0;
            float proj_this_min_t = merged_edge.proj_a_point(m_inter_point + m_end_min * m_inter_direction);
            Eigen::Vector3f pt_test1 = m_inter_point + m_end_min * m_inter_direction,
                            pt_test2 = merged_edge.m_inter_point + proj_this_min_t * merged_edge.m_inter_direction;
            float ratio_test = proj_a_point(pt_test2);
            // std::cout << "before proj this min: " << pt_test1.transpose() << ", after proj: " << pt_test2.transpose() 
            //           << ", ratio before proj=" << m_end_min << ", ratio after proj=" << ratio_test<< "\n";
            
            float proj_this_max_t = merged_edge.proj_a_point(m_inter_point + m_end_max * m_inter_direction);
            // pt_test1 = m_inter_point + m_end_max * m_inter_direction;
            // pt_test2 = merged_edge.m_inter_point + proj_this_max_t * merged_edge.m_inter_direction;
            // ratio_test = proj_a_point(pt_test2);
            // std::cout << "before proj this max: " << pt_test1.transpose() << ", after proj: " << pt_test2.transpose() 
            //             << ", ratio before proj=" << m_end_max << ", ratio after proj=" << ratio_test<< "\n";
            
            float proj_other_min_t = merged_edge.proj_a_point(other.m_inter_point + other.m_end_min * other.m_inter_direction);
            // pt_test1 = other.m_inter_point + other.m_end_min * other.m_inter_direction; 
            // pt_test2 = merged_edge.m_inter_point + proj_other_min_t * merged_edge.m_inter_direction;
            // ratio_test = other.proj_a_point(pt_test2);
            // std::cout << "before proj other min: " << pt_test1.transpose() << ", after proj: " << pt_test2.transpose() 
            //             << ", ratio before proj=" << other.m_end_min << ", ratio after proj=" << ratio_test<< "\n";
            
            float proj_other_max_t = merged_edge.proj_a_point(other.m_inter_point + other.m_end_max * other.m_inter_direction);
            // pt_test1 = other.m_inter_point + other.m_end_max * other.m_inter_direction;
            // pt_test2 = merged_edge.m_inter_point + proj_other_max_t * merged_edge.m_inter_direction;
            // std::cout << "before proj other max: " << pt_test1.transpose() << ", after proj: " << pt_test2.transpose() 
            //             << ", ratio before proj=" << other.m_end_max << ", ratio after proj=" << ratio_test<< "\n";

            float this_max_t = std::max(proj_this_min_t, proj_this_max_t);
            float this_min_t = std::min(proj_this_min_t, proj_this_max_t);
            float other_max_t = std::max(proj_other_min_t, proj_other_max_t);
            float other_min_t = std::min(proj_other_min_t, proj_other_max_t);

            float minmax_end = std::min(this_max_t, other_max_t);
            float maxmin_end = std::max(this_min_t, other_min_t);
            float overlap = minmax_end - maxmin_end;
            if (overlap <= dist_thresh) return 0;

            merged_edge.m_end_min = maxmin_end;
            merged_edge.m_end_max = minmax_end;
            merged_edge.update_two_ends(maxmin_end, minmax_end);

            new_edges.push_back(merged_edge);

            if (this_max_t - minmax_end > dist_thresh || maxmin_end - this_min_t > dist_thresh)
            {
                // TODO: DE THE BUG HERE
                // project back to the original edges
                // preserve the remaining original edges: could be 0, 1, 2 edges
                float proj_overlap_max_on_this = proj_a_point(merged_edge.m_inter_point + minmax_end * merged_edge.m_inter_direction);
                proj_overlap_max_on_this = std::clamp(proj_overlap_max_on_this, m_end_min, m_end_max);
                float proj_overlap_min_on_this = proj_a_point(merged_edge.m_inter_point + maxmin_end * merged_edge.m_inter_direction);
                proj_overlap_min_on_this = std::clamp(proj_overlap_min_on_this, m_end_min, m_end_max);

                float overlap_max_on_this = std::max(proj_overlap_max_on_this, proj_overlap_min_on_this);
                float overlap_min_on_this = std::min(proj_overlap_max_on_this, proj_overlap_min_on_this);
                
                if (m_end_max - overlap_max_on_this > dist_thresh)
                {
                    Edge remain_edge = Edge({m_plane1, m_plane2}, {m_flip1, m_flip2}, m_inward_radian,
                                            m_inter_direction, m_inter_point, m_proj_from_origin, m_distance_from_origin,
                                            overlap_max_on_this, m_end_max);
                    new_edges.push_back(remain_edge);
                }
                if (overlap_min_on_this - m_end_min > dist_thresh)
                {
                    Edge remain_edge = Edge({m_plane1, m_plane2}, {m_flip1, m_flip2}, m_inward_radian,
                                            m_inter_direction, m_inter_point, m_proj_from_origin, m_distance_from_origin,
                                            m_end_min, overlap_min_on_this);
                    new_edges.push_back(remain_edge);
                }
            }

            if (other_max_t - minmax_end > dist_thresh || maxmin_end - other_min_t > dist_thresh)
            {
                // TODO: DE THE BUG HERE
                float proj_overlap_max_on_other = other.proj_a_point(merged_edge.m_inter_point + minmax_end * merged_edge.m_inter_direction);
                proj_overlap_max_on_other = std::clamp(proj_overlap_max_on_other, other.m_end_min, other.m_end_max);
                float proj_overlap_min_on_other = other.proj_a_point(merged_edge.m_inter_point + maxmin_end * merged_edge.m_inter_direction);
                proj_overlap_min_on_other = std::clamp(proj_overlap_min_on_other, other.m_end_min, other.m_end_max);
                float overlap_max_on_other = std::max(proj_overlap_max_on_other, proj_overlap_min_on_other);
                float overlap_min_on_other = std::min(proj_overlap_max_on_other, proj_overlap_min_on_other);
                
                if (other.m_end_max - overlap_max_on_other > dist_thresh)
                {
                    Edge remain_edge = Edge({other.m_plane1, other.m_plane2}, {other.m_flip1, other.m_flip2}, other.m_inward_radian,
                                            other.m_inter_direction, other.m_inter_point, other.m_proj_from_origin, other.m_distance_from_origin,
                                            overlap_max_on_other, other.m_end_max);
                    new_edges.push_back(remain_edge);
                }
                if (overlap_min_on_other - other.m_end_min > dist_thresh)
                {
                    Edge remain_edge = Edge({other.m_plane1, other.m_plane2}, {other.m_flip1, other.m_flip2}, other.m_inward_radian,
                                            other.m_inter_direction, other.m_inter_point, other.m_proj_from_origin, other.m_distance_from_origin,
                                            other.m_end_min, overlap_min_on_other);
                    new_edges.push_back(remain_edge);
                }
            }
            return 1;
        }
        return 0;
    }

    float proj_a_point(const Eigen::Vector3f &pt) const
    {
        return (pt.dot(m_inter_direction) - m_inter_point.dot(m_inter_direction)) / m_inter_direction.dot(m_inter_direction);
    }

    float proj_a_point(const Point &pt) const
    {
        Eigen::Vector3f v(
            CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()));
        return proj_a_point(v);
    }

    void update_parameter_of_intersection()
    {
        float a1 = CGAL::to_double(m_plane1->a()), b1 = CGAL::to_double(m_plane1->b()), c1 = CGAL::to_double(m_plane1->c());
        float a2 = CGAL::to_double(m_plane2->a()), b2 = CGAL::to_double(m_plane2->b()), c2 = CGAL::to_double(m_plane2->c());
        Eigen::Vector3f n1(a1, b1, c1), n2(a2, b2, c2);
        float d1 = CGAL::to_double(m_plane1->d()), d2 = CGAL::to_double(m_plane2->d());
        m_inter_direction = n1.cross(n2);
        m_inter_direction.normalize();
        float ai = m_inter_direction(0), bi = m_inter_direction(1), ci = m_inter_direction(2);
       
        Eigen::Matrix2f A;
        Eigen::Vector2f b = {-d1, -d2};
        if (std::abs(ai) >= std::max(std::abs(bi), std::abs(ci)))
        {
            A << b1, c1, b2, c2;
            Eigen::Vector2f x = A.inverse() * b;
            m_inter_point = {0, x(0), x(1)};
        }
        else if (std::abs(bi) >= std::max(std::abs(ci), std::abs(ai)))
        {
            A << c1, a1, c2, a2;
            Eigen::Vector2f x = A.inverse() * b;
            m_inter_point = {x(1), 0, x(0)};
        }
        else // (std::abs(ci) > std::max(std::abs(ai), std::abs(bi)))
        {
            A << a1, b1, a2, b2;
            Eigen::Vector2f x = A.inverse() * b;
            m_inter_point = {x(0), x(1), 0};
        }
        
        // check nan
        if (m_inter_direction(0) != m_inter_direction(0) || m_inter_direction(1) != m_inter_direction(1) || m_inter_direction(2) != m_inter_direction(2)
         || m_inter_point(0) != m_inter_point(0) || m_inter_point(1) != m_inter_point(1) || m_inter_point(2) != m_inter_point(2))
            std::cout << "m_inter_direction=" << m_inter_direction.transpose() << ", m_inter_point=" << m_inter_point.transpose() << "\n";
    }

    std::vector<Edge> split(const LCC_3 &lcc) const
    {
        std::vector<Edge> split_edges;
        std::map<Edge_idx, bool> is_join;
        for (auto [eid, e_dart] : m_cm_edges)
            is_join[eid] = false;

        for (auto [eid, e_dart] : m_cm_edges)
        {
            if (is_join[eid]) continue;
            is_join[eid] = true;
            Vertex_idx end1 = lcc.info<0>(e_dart), end2 = lcc.info<0>(lcc.beta<1>(e_dart));
            Dart_const_descriptor end1_dart = e_dart, end2_dart = lcc.beta<1>(e_dart);

            bool close_end1 = false, close_end2 = false;
            do {
                close_end1 = true;
                close_end2 = true;
                for (auto [eid2, e_dart2] : m_cm_edges)
                {
                    if (is_join[eid2]) continue;
                    Vertex_idx end3 = lcc.info<0>(e_dart2), end4 = lcc.info<0>(lcc.beta<1>(e_dart2));
                    if (end3 == end1 || end4 == end1)
                    {
                        is_join[eid2] = true;
                        close_end1 = false;
                        end1 = end3 == end1 ? end4 : end3;
                        end1_dart = end3 == end1 ? lcc.beta<1>(e_dart2) : e_dart2;
                        break;
                    }
                    else if (end3 == end2 || end4 == end2)
                    {
                        is_join[eid2] = true;
                        close_end2 = false;
                        end2 = end3 == end2 ? end4 : end3;
                        end2_dart = end3 == end2 ? lcc.beta<1>(e_dart2) : e_dart2;
                        break;
                    }
                }

            }while(!close_end1 && !close_end2);

            Edge new_edge = Edge({m_plane1, m_plane2}, {m_flip1, m_flip2}, m_inward_radian,
                                  m_inter_direction, m_inter_point, m_proj_from_origin, m_distance_from_origin,
                                  lcc.point(end1_dart), lcc.point(end2_dart));
            split_edges.push_back(new_edge);
        }

        return split_edges;
    }

    void update_two_ends(const LCC_3 &lcc)
    {
        all_t.clear();
        all_ends.clear();
        for (auto [eid, e_dart] : m_cm_edges)
        {
            Point p1 = lcc.point(e_dart), p2 = lcc.point(lcc.beta<1>(e_dart));
            all_ends.push_back(p1);
            all_ends.push_back(p2);
            float t1 = proj_a_point(p1);
            float t2 = proj_a_point(p2);
            all_t.push_back(t1);
            all_t.push_back((t1 + t2) / 2);
            all_t.push_back(t2);
            Eigen::Vector3f mid_pt = m_inter_point + ((t1+t2)/2) * m_inter_direction;
            // std::cout << "m_inter_point=" << m_inter_point.transpose()
            //           << "m_inter_direction=" << m_inter_direction.transpose() 
            //           << "mid_pt=" << mid_pt.transpose() << "\n";
            all_ends.push_back(Point(mid_pt(0), mid_pt(1), mid_pt(2)));
        }
        m_end_max = *std::max_element(all_t.begin(), all_t.end());
        m_end_min = *std::min_element(all_t.begin(), all_t.end());
        m_length = m_end_max - m_end_min;
    }

    void update_two_ends(
        const Point &end1, const Point &end2)
    {
        Eigen::Vector3f v1(
            CGAL::to_double(end1.x()), CGAL::to_double(end1.y()), CGAL::to_double(end1.z()));
        Eigen::Vector3f v2(
            CGAL::to_double(end2.x()), CGAL::to_double(end2.y()), CGAL::to_double(end2.z()));
        update_two_ends(v1, v2);
    }

    void update_two_ends(
        const Eigen::Vector3f end1, const Eigen::Vector3f end2)
    {
        float t1 = (end1.dot(m_inter_direction) - m_inter_point.dot(m_inter_direction))
                        / (m_inter_direction.dot(m_inter_direction));
        float t2 = (end2.dot(m_inter_direction) - m_inter_point.dot(m_inter_direction))
                    / (m_inter_direction.dot(m_inter_direction));
        
        // put two ends and the middle point in the all_t and all_ends
        all_t.clear();
        
        all_t.push_back(t1);
        all_t.push_back((t1 + t2) / 2);
        all_t.push_back(t2);

        all_ends.clear();
        Point p1(end1(0), end1(1), end1(2)), p2(end2(0), end2(1), end2(2));
        Point mid((end1(0) + end2(0)) / 2, (end1(1) + end2(1)) / 2, (end1(2) + end2(2)) / 2);
        all_ends.push_back(p1);
        all_ends.push_back(mid);
        all_ends.push_back(p2);

        m_end_min = *std::min_element(all_t.begin(), all_t.end());
        m_end_max = *std::max_element(all_t.begin(), all_t.end());
        m_length = m_end_max - m_end_min;
    }

    void update_two_ends(
        const float end_min, const float end_max)
    {
        if (end_min >= end_max) 
            std::cout << "end_min=" << end_min << ", end_max=" << end_max << "\n";
        assert(end_min < end_max);
        m_end_min = end_min;
        m_end_max = end_max;
        m_length = m_end_max - m_end_min;

        all_t.clear();
        all_t.push_back(m_end_min);
        all_t.push_back((m_end_min + m_end_max) / 2);
        all_t.push_back(m_end_max);

        all_ends.clear();
        Eigen::Vector3f v1 = m_inter_point + m_end_min * m_inter_direction,
                        v2 = m_inter_point + m_end_max * m_inter_direction,
                        v3 = m_inter_point + (m_end_min + m_end_max) / 2 * m_inter_direction;
        Point p1 = Point(v1(0), v1(1), v1(2)),
              p2 = Point(v2(0), v2(1), v2(2)),
              p3 = Point(v3(0), v3(1), v3(2));
        all_ends.push_back(p1);
        all_ends.push_back(p2);
        all_ends.push_back(p3);
    }

    void update_proj_from_origin()
    {
        m_proj_from_origin = - (m_inter_point.dot(m_inter_direction)) / m_inter_direction.dot(m_inter_direction);
        Eigen::Vector3f proj_pt = m_inter_point + m_proj_from_origin * m_inter_direction;
        // std::cout << "proj_pt.dot(inter_direction)=" << proj_pt.dot(m_inter_direction) << "\n";
        m_distance_from_origin = proj_pt.norm();
    }

    void update_weight(const Tree &ch_tree)
    {
        distances_from_ends_to_ch.resize(all_ends.size(), 0);
        Eigen::Vector3f n1 (
            CGAL::to_double(m_plane1->a()),
            CGAL::to_double(m_plane1->b()),
            CGAL::to_double(m_plane1->c()));
        Eigen::Vector3f n2 (
            CGAL::to_double(m_plane2->a()),
            CGAL::to_double(m_plane2->b()),
            CGAL::to_double(m_plane2->c()));
        if (m_flip1) n1 = -n1;
        if (m_flip2) n2 = -n2;
        Eigen::Vector3f out_bisector_dir = (n1 + n2); // should be -(n1+n2), but that doesn't work
        Vector ray_dir(out_bisector_dir(0), out_bisector_dir(1), out_bisector_dir(2));
        float total_dist = 0;
        int count_dist = 0;
        for(size_t i = 0; i < all_ends.size(); i++)
        {
            Ray ray(all_ends[i], ray_dir);
            std::list<Ray_intersection> intersections;
            ch_tree.all_intersections(ray, std::back_inserter(intersections));
            bool intersected_a_point = false;
            for (auto &inter : intersections)
            {
                if (inter)
                {
                    if (std::get_if<Point>(&(inter->first)))
                    {
                        const Point* p =  std::get_if<Point>(&(inter->first));
                        float dist = std::sqrt(CGAL::to_double(CGAL::squared_distance(*p, all_ends[i])));
                        total_dist += dist;
                        count_dist++;
                        intersected_a_point = true;
                        break;
                    }
                }
            }
            // if (!intersected_a_point)
            //     std::cerr << "no intersection point found\n";
        }
        if (count_dist > 0)
            m_distance_to_ch = total_dist / count_dist;
        else
            m_distance_to_ch = 0;
        m_weight = m_distance_to_ch * m_length;
    }

    const Plane* m_plane1;
    const Plane* m_plane2;
    bool m_flip1;
    bool m_flip2;
    float m_inward_radian;
    // parameters of the intersection line
    Eigen::Vector3f m_inter_direction;
    Eigen::Vector3f m_inter_point; 
    std::vector<Point> all_ends;
    std::vector<float> all_t;
    float m_length;
    float m_end_min;
    float m_end_max;
    float m_proj_from_origin;
    float m_distance_from_origin;
    float m_weight;
    std::vector<float> distances_from_ends_to_ch;
    float m_distance_to_ch;
    // set of the cm edges on the intersection line
    std::map<Edge_idx, Dart_const_descriptor> m_cm_edges;
};

class EdgeSet{

public: 
    EdgeSet(){}

    size_t size() const
    {
        return m_edges.size();
    }

    void clear() 
    {
        m_edges.clear();
    }

    void insert(
        const std::pair<const Plane*, const Plane*> &edge,
        const std::pair<bool, bool> &flip,
        const Edge_idx edge_idx,
        Dart_const_descriptor edge_dart,
        const float inward_radian)
    {
        bool inserted = false;
        for (auto &me : m_edges)
        {
            if (edge.first == me.m_plane1 && edge.second == me.m_plane2
               && flip.first == me.m_flip1 && flip.second == me.m_flip2)
            {
                me.insert_cm_edge(edge_idx, edge_dart);
                inserted = true;
                break;
            }
        }
        if (!inserted)
        {
            Edge new_edge(edge, flip, inward_radian);
            if (new_edge.m_inter_direction.norm() > 0) 
            {
                new_edge.insert_cm_edge(edge_idx, edge_dart);
                m_edges.push_back(new_edge);
            }
        }
    }

    void merge_edges(const LCC_3 &lcc, Volume_idx cell_idx, const float dist_thresh)
    {
        bool updated = false;
        std::cout << "begin to merge...\n";
        do {
            updated = false;
            for (auto it1 = m_edges.begin(); it1 != m_edges.end(); ++it1)
            {
                // std::cout << "it1=" << std::distance(m_edges.begin(), it1) << '\n';
                for (auto it2 = std::next(it1); it2 != m_edges.end(); ++it2)
                {
                    std::vector<Edge> new_edges;
                    int result = it1->merge(lcc, cell_idx, *it2, dist_thresh, new_edges);
                    if (result == 0) continue;
                    else if (result == 1)
                    {
                        assert(new_edges.size() > 0);
                        m_edges.erase(it2);
                        // std::cout << "after erasing it1=" << std::distance(m_edges.begin(), it1) << '\n'; 
                        m_edges.erase(it1);
                        m_edges.insert(m_edges.end(), new_edges.begin(), new_edges.end());
                        updated = true;
                        break;
                    }
                    else
                    {
                        std::cout << "result==2\n";
                        m_edges.erase(it2);
                        m_edges.erase(it1);
                        updated = true;
                        break;
                    }
                }
                if (updated) break;
            }
        }while(updated);
        std::cout << "merge done!\n";
    }

    void update_edge_ends(const LCC_3 &lcc)
    {
        for (auto &edge: m_edges)
            edge.update_two_ends(lcc);
    }

    void split_edges(const LCC_3 &lcc)
    {
        std::cout << "#edge before splitting=" << m_edges.size();

        std::vector<Edge> updated_edges;
        for (auto &edge : m_edges)
        {
            std::vector<Edge> split_edges = edge.split(lcc);
            updated_edges.insert(updated_edges.end(), split_edges.begin(), split_edges.end());
        }
        m_edges = updated_edges;
        std::cout << ", after splitting=" << m_edges.size() << '\n';
    }

    void update_weights(const Tree &ch_tree)
    {
        m_max_weight = 0;
        m_min_weight = std::numeric_limits<float>::max();
        for (auto &edge : m_edges)
        {
            edge.update_weight(ch_tree);
            // std::cout << "updated weight=" << edge.m_weight 
            //           << "(" << edge.m_distance_to_ch << "*" << edge.m_length << ")\n";
            m_max_weight = std::max(m_max_weight, edge.m_weight);
            m_min_weight = std::min(m_min_weight, edge.m_weight);
        }

        m_edges.erase(
            std::remove_if(m_edges.begin(), m_edges.end(),
                [this](const Edge &edge) -> bool
                {
                    return edge.m_weight < 1e-6;
                }), m_edges.end()
        );
    }

    void filter_edges(const float lod)
    {
        float dx = CGAL::to_double(Global_var::box.xmax() - Global_var::box.xmin());
        float dy = CGAL::to_double(Global_var::box.ymax() - Global_var::box.ymin());
        float dz = CGAL::to_double(Global_var::box.zmax() - Global_var::box.zmin());

        float vertical_diagonal_face_area = dz * std::sqrt(dx * dx + dy * dy);
        
        m_edges.erase(
            std::remove_if(m_edges.begin(), m_edges.end(),
                [lod, vertical_diagonal_face_area](const Edge &edge) -> bool
                {
                    return edge.m_weight < lod * vertical_diagonal_face_area;
                }), m_edges.end()
        );
    }

    void remove_convex_edges(
        const float concave_radian_thresh)
    {
        m_edges.erase(
            std::remove_if(m_edges.begin(), m_edges.end(),
                [concave_radian_thresh](const Edge &edge) -> bool
                {
                    return edge.m_inward_radian < concave_radian_thresh;
                }), m_edges.end());
    }

    void sort()
    {
        std::sort(m_edges.begin(), m_edges.end(),
            [this](const Edge &e1, const Edge&e2)
            {
                return e1.m_weight > e2.m_weight;
            });
    }

    std::vector<Edge> m_edges;
    float m_max_weight;
    float m_min_weight;
};