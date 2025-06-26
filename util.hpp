#pragma once

#include <algorithm>
#include <filesystem>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <ctime>
#include <list>
#include <memory>
#define CGAL_USE_BASIC_VIEWER

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
// #include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/algorithm.h>
#include <CGAL/bounding_box.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>


#if VIS == 1
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/draw_surface_mesh.h>
#include<CGAL/draw_arrangement_2.h>
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

// rapidjson includes
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>

#include <sdf/sdf.hpp>

namespace rj = rapidjson;

std::array<std::array<int, 3>, 50> pastelColorMap = {{
    {255, 0, 0},     // Red
    {0, 255, 0},     // Green
    {0, 0, 255},     // Blue
    {255, 255, 0},   // Yellow
    {255, 0, 255},   // Magenta
    {0, 255, 255},   // Cyan
    {255, 128, 0},   // Orange
    {128, 0, 255},   // Purple
    {255, 128, 128}, // Light Pink
    {128, 255, 128}, // Light Green
    {128, 128, 255}, // Light Blue
    {255, 255, 128}, // Light Yellow
    {255, 128, 255}, // Light Magenta
    {128, 255, 255}, // Light Cyan
    {192, 0, 0},     // Dark Red
    {0, 192, 0},     // Dark Green
    {0, 0, 192},     // Dark Blue
    {192, 192, 0},   // Dark Yellow
    {192, 0, 192},   // Dark Magenta
    {0, 192, 192},   // Dark Cyan
    {255, 64, 64},   // Bright Red
    {64, 255, 64},   // Bright Green
    {64, 64, 255},   // Bright Blue
    {255, 255, 64},  // Bright Yellow
    {255, 64, 255},  // Bright Magenta
    {64, 255, 255},  // Bright Cyan
    {128, 128, 0},   // Olive
    {128, 0, 128},   // Indigo
    {0, 128, 128},   // Teal
    {255, 128, 64},  // Bright Orange
    {128, 64, 255},  // Bright Purple
    {64, 255, 128},  // Bright Lime
    {128, 255, 64},  // Bright Chartreuse
    {64, 128, 255},  // Bright Azure
    {255, 64, 128},  // Bright Fuchsia
    {128, 64, 0},    // Brown
    {0, 128, 64},    // Forest Green
    {64, 0, 128},    // Navy
    {255, 165, 0},   // Gold
    {173, 216, 230}, // LightBlue
    {240, 128, 128}, // LightCoral
    {50, 205, 50},   // LimeGreen
    {123, 104, 238}, // MediumSlateBlue
    {255, 20, 147},  // DeepPink
    {72, 61, 139},   // DarkSlateBlue
    {176, 196, 222}, // LightSteelBlue1
    {255, 228, 196}, // Bisque
    {245, 222, 179}, // Wheat
    {255, 240, 245}, // LavenderBlush
    {224, 255, 255}, // LightCyan
}};

std::array<int, 3> hsv_to_rgb(float h, float s, float v) {
    float r, g, b;
    int i = static_cast<int>(h / 60.0f) % 6;
    float f = h / 60.0f - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);

    switch (i) {
    case 0:
        r = v;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = v;
        b = p;
        break;
    case 2:
        r = p;
        g = v;
        b = t;
        break;
    case 3:
        r = p;
        g = q;
        b = v;
        break;
    case 4:
        r = t;
        g = p;
        b = v;
        break;
    case 5:
        r = v;
        g = p;
        b = q;
        break;
    }

    return {static_cast<int>(r * 255), static_cast<int>(g * 255),
            static_cast<int>(b * 255)};
}

void create_colormap(int num_colors,
            std::vector<CGAL::IO::Color> &colors) {
    colors.clear();
    std::vector<std::tuple<int, int, int>> rgb_colors;
    for (int i = 0; i < num_colors; ++i) {
        float h = i * (360.0f / float(num_colors)); // 色相从0到360度
        float s = 1.0f;                             // 饱和度为100%
        float v = 1.0f;                             // 亮度为100%
        std::array<int, 3> rgb = hsv_to_rgb(h, s, v);
        colors.push_back(CGAL::IO::Color(rgb[0], rgb[1], rgb[2]));
    }
}

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
// typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
typedef Kernel::FT FT;
typedef Kernel::RT RT;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_2 Point2;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Vector_2 Vector2;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Ray_3 Ray;
typedef Kernel::Triangle_3 Triangle;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> Traits;

typedef Surface_mesh::Face_index Face_index;
typedef Surface_mesh::Vertex_index Vertex_index;
typedef Surface_mesh::Edge_index Edge_index;
typedef Surface_mesh::Halfedge_index Halfedge_index;
typedef Surface_mesh::Point Mesh_point;
typedef CGAL::Partition_traits_2<Kernel, CGAL::Pointer_property_map<Point2>::type > Partition_traits_2; 
typedef Partition_traits_2::Polygon_2                       Polygon2;  // a polygon of indices

using Vertex_idx = int;
using Edge_idx = int;
using Face_idx = int;
using Volume_idx = int;
namespace PMP = CGAL::Polygon_mesh_processing;

struct Items_with_indices {
    // using Use_index=CGAL::Tag_true; // use indices
    // using Index_type=std::uint16_t; // 16 bits
    template <class Refs> struct Dart_wrapper {
        typedef CGAL::Cell_attribute<Refs, std::pair<int, bool>> Vol_attrib; // explicit indexing
        typedef CGAL::Cell_attribute<Refs, std::pair<int, const Plane*>>
            Face_attrib; // explicit indexing
        typedef CGAL::Cell_attribute<Refs, int>
            Edge_attrib; // explicit indexing
        typedef CGAL::Cell_attribute_with_point<Refs, int, CGAL::Tag_true>
            Vtx_attrib; // explicit indexing
        typedef std::tuple<Vtx_attrib, Edge_attrib, Face_attrib, Vol_attrib>
            Attributes;
    };
};

// typedef CGAL::Linear_cell_complex_for_generalized_map<3, 3, Traits,
//                                                       Items_with_indices>
typedef CGAL::Linear_cell_complex_for_combinatorial_map<
    3, 3, Traits, Items_with_indices> LCC_3;
typedef LCC_3::Dart_descriptor Dart_descriptor;
typedef LCC_3::Dart_const_descriptor Dart_const_descriptor;
typedef LCC_3::Attribute_type<3>::type Vol_attribute;
typedef LCC_3::Attribute_type<2>::type Face_attribute;
typedef LCC_3::Attribute_type<1>::type Edge_attribute;
typedef LCC_3::Attribute_type<0>::type Vtx_attribute;
typedef LCC_3::size_type size_type;

struct Split_functor_3 {
    Split_functor_3(LCC_3 &lcc) : lcc_(lcc) {}
    void operator()(Vol_attribute &ca1, Vol_attribute &ca2) {
        // lcc_.set_attribute<3>(ca2.dart(), lcc_.create_attribute<3>()); //
        // create new attribute automatically?
        // int top_3 = -1;
        // for (typename LCC_3::Dart_range::iterator itall = lcc_.darts().begin(),
        //                                           itallend = lcc_.darts().end();
        //      itall != itallend; ++itall) {
        //     top_3 = std::max(top_3, lcc_.info<3>(itall));
        // }
        lcc_.info<3>(ca2.dart()).first = lcc_.number_of_attributes<3>()-1;
        lcc_.info<3>(ca2.dart()).second = lcc_.info<3>(ca1.dart()).second;
    }

  private:
    LCC_3 &lcc_;
};

struct Split_functor_2 {
    Split_functor_2(LCC_3 &lcc) : lcc_(lcc) {}
    void operator()(Face_attribute &ca1, Face_attribute &ca2) {
        // assume face splitting is always on the same plane
        // assign same plane to the split face
        lcc_.info<2>(ca2.dart()) = std::make_pair(
                lcc_.number_of_attributes<2>() - 1,
                lcc_.info<2>(ca1.dart()).second); 
    }

  private:
    LCC_3 &lcc_;
};

struct Split_functor_1 {
    Split_functor_1(LCC_3 &lcc) : lcc_(lcc) {}
    void operator()(Edge_attribute &ca1, Edge_attribute &ca2) {
        // lcc_.set_attribute<1>(ca2.dart(), lcc_.create_attribute<1>());
        lcc_.info<1>(ca2.dart()) = lcc_.number_of_attributes<1>() - 1;
        // int top_1 = -1;
        // for (typename LCC_3::Dart_range::iterator itall =
        // lcc_.darts().begin(),
        //    itallend = lcc_.darts().end(); itall!=itallend; ++itall)
        // {
        //     top_1 = std::max(top_1, lcc_.info<3>(itall));
        // }
        // lcc_.info<1>(ca2.dart()) = top_1 + 1;
    }

  private:
    LCC_3 &lcc_;
};

// Functor called when two faces are merged.
struct Merge_functor_3
{
  // operator() automatically called before a merge.
  void operator()(Vol_attribute& ca1, Vol_attribute& ca2)
  {
    // std::cout   <<"Before on merge volumes: " 
    //             << "volume1=" << ca1.info().first << "(in: " << ca1.info().second << "), " 
    //             << " volume2="<<ca2.info().first << "(in: " << ca2.info().second << ")\n";

    ca1.info().first=std::min(ca1.info().first, ca2.info().first); // Update can be done incrementally.
    assert (ca1.info().second == ca2.info().second);
  }
};

struct Merge_functor_2
{
  // operator() automatically called before a merge.
  void operator()(Face_attribute &ca1, Face_attribute &ca2)
  {
    // std::cout   <<"Before on merge volumes: " 
    //             << "volume1=" << ca1.info().first << "(in: " << ca1.info().second << "), " 
    //             << " volume2="<<ca2.info().first << "(in: " << ca2.info().second << ")\n";

    ca1.info().first=std::min(ca1.info().first, ca2.info().first); // Update can be done incrementally.
    assert (ca1.info().second == ca2.info().second);
  }
};

// struct Draw_functor : public CGAL::draw_function_for_lcc
// {
//     Draw_functor(LCC_3::size_type fm) :
//             is_highlighted(fm), highlight_f(true) {}

//     template<typename LCC>
//     bool colored_face(const LCC& /* alcc */,
//                     typename LCC::Dart_const_descriptor /* d */) const {if
//                     (highlight_f) return true; else return false;}

//     template<typename LCC>
//     CGAL::IO::Color face_color(const LCC&  alcc ,
//                             typename LCC::Dart_const_descriptor d) const
//     {
//         if (highlight_f)
//         {
//             if (alcc.is_marked(d, is_highlighted)) return
//             CGAL::IO::Color(255, 0, 0); // in red else return
//             CGAL::IO::Color(0, 255, 0); // in green
//         }
//     }

//     LCC_3::size_type is_highlighted;
//     bool highlight_f;
// };

#if VIS == 1
struct LCC_gs_options :
  public CGAL::Graphics_scene_options<LCC_3,
                               typename LCC_3::Dart_const_descriptor,
                               typename LCC_3::Dart_const_descriptor,
                               typename LCC_3::Dart_const_descriptor,
                               typename LCC_3::Dart_const_descriptor>
{
  LCC_gs_options(LCC_3::size_type highlight) : m_highlight(highlight){}

  bool colored_face(const LCC_3& lcc, const typename LCC_3::Dart_const_descriptor d) const
  {
    return true;
  }

  CGAL::IO::Color face_color(const LCC_3& lcc, const typename LCC_3::Dart_const_descriptor d) const
  {
    if (lcc.is_marked(d, m_highlight)) 
        return CGAL::IO::Color(255, 0, 0);

    return CGAL::IO::Color(0, 255, 0);
  }

  LCC_3::size_type m_highlight;
};
#endif

// struct Draw_face_functor : public CGAL::DefaultDrawingFunctorLCC
// {
//     Draw_face_functor(LCC_3::size_type fm, LCC_3::size_type h) :
//     show_face(fm), is_highlighted(fm), highlight_f(true) {}
//     Draw_face_functor(LCC_3::size_type vm) : show_face(vm),
//     highlight_f(false) {}

//     template<typename LCC>
//     bool colored_face(const LCC& /* alcc */,
//                     typename LCC::Dart_const_descriptor /* d */) const {if
//                     (highlight_f) return true; else return false;}

//     template<typename LCC>
//     CGAL::IO::Color face_color(const LCC&  alcc ,
//                             typename LCC::Dart_const_descriptor d) const
//     {
//         if (highlight_f)
//         {
//             if (alcc.is_marked(d, is_highlighted)) return
//             CGAL::IO::Color(255, 0, 0); // in green else return
//             CGAL::IO::Color(0, 255, 0); // out red
//         }
//     }

//     template<typename LCC>
//     bool draw_face(const LCC& alcc,
//                    typename LCC::Dart_const_descriptor d) const
//     {
//         if (alcc.is_marked(d, show_face)) return true;
//         else return false;
//     }

//     LCC_3::size_type show_face;
//     LCC_3::size_type is_highlighted;
//     bool highlight_f;
// };

using Neighbor_query =
    CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Surface_mesh>;
using Region_type =
    CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<
        Kernel, Surface_mesh>;
using Sorting =
    CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<
        Kernel, Surface_mesh, Neighbor_query>;
using Region_growing =
    CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;
using Region = std::vector<Face_index>;
using Primitive_and_region = std::pair<typename Region_type::Primitive, Region>;

// using Item = typename RegionType::Item;
using Region = std::vector<Face_index>;
using Primitive_and_region = std::pair<typename Region_type::Primitive, Region>;

typedef Kernel::Iso_cuboid_3 BBox;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits_3<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;
typedef std::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;


typedef CGAL::AABB_triangle_primitive_3<Kernel, std::list<Triangle>::const_iterator> Primitive_tri;
typedef CGAL::AABB_traits_3<Kernel, Primitive_tri> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree_tri;
typedef std::optional< Tree_tri::Intersection_and_primitive_id<Plane>::Type > Plane_intersection;

typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef Segment_traits_2::X_monotone_curve_2              Segment_2;
typedef CGAL::Arr_consolidated_curve_data_traits_2<Segment_traits_2, Dart_descriptor>  Traits_2;
typedef Traits_2::X_monotone_curve_2                      CM_Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                     Arrangement_2;

typedef CGAL::Triangulation_3<Kernel>      Triangulation;
typedef Triangulation::Cell_handle     Cell_handle;
typedef Triangulation::Vertex_handle   Vertex_handle;
typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Finite_cells_iterator    Finite_cells_iterator;

using In_cell_pts = std::vector<int>;

using SDF_Points = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
using SDF_Triangles = Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using SDF_Result = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

class Component;
using Evaluator =
    std::function<float(const Component *, const std::vector<Component *> &)>;

struct Concave_comp {
    // Operator() overloading
    bool operator()(const std::pair<const Plane *, const Plane *> &e1,
                    const std::pair<const Plane *, const Plane *> &e2) const {
        // increasing order, the larger address, the larger index in the set
        if (e1.first > e2.first ||
            (e1.first == e2.first && e1.second > e2.second))
            return true;
        else
            return false;
    }
};

enum class Cut_mode { A, H, V, HV, HI, VI};

float volume_of_tetrahedron(
    const std::vector<Eigen::Vector3f> &four_pts)
{
    float volume = 0;
    for (size_t i = 0 ; i < 4; i++)
    {
        Eigen::Vector3f p0 = four_pts[i],
                        p1 = four_pts[(i+1)%4],
                        p2 = four_pts[(i+2)%4];
        volume += 1.0/6.0 * 
                (-p2(0) * p1(1) * p0(2) + p1(0) * p2(1) * p0(2) + p2(0) * p0(1) * p1(2) 
                - p0(0) * p2(1) * p1(2) - p1(0) * p0(1) * p2(2) + p0(0) * p1(1) * p2(2));
    }

    return std::abs(volume);
}

Eigen::MatrixXf make_para_comb(
    const int num_para,
    const int num_pins)
{
    std::vector<float> options (num_pins, 0);
    for (size_t i = 0; i < num_pins; i++)
        options[i] = (i+1) *  1.0 / (num_pins + 1);
    
    size_t num_combs = std::pow(num_pins, num_para-1);
    // given num_param = 4, -> lambda1, lambda2, lambda3, lambda4 = 1 - lambda1 - lambda2 - lambda3
    Eigen::MatrixXf combs(num_combs, num_para);
    for (int i = 0; i < num_combs; i++)
    {
        float accu = 0;
        for (int j = 0; j < num_para-1; j++)
        {
            combs(i, j) = (1 - accu) * options[int(i / std::pow(num_pins, num_para-2-j)) % num_pins];
            accu += combs(i, j);
        }
        combs(i, num_para-1) = 1 - accu;
    }
    std::cout << "combs: " << combs << '\n';
    return combs;
}

BBox bbox_of_cells_3(const std::map<int, Dart_descriptor> &cells, LCC_3 &lcc) {
    std::vector<Point> pts;
    for (const auto &[idx, dart] : cells) {
        for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator v =
                 lcc.one_dart_per_incident_cell<0, 3>(dart).begin();
             v != lcc.one_dart_per_incident_cell<0, 3>(dart).end(); ++v)
            pts.push_back(lcc.point(v));
    }
    return CGAL::bounding_box(pts.begin(), pts.end());
}

Vector move_mesh_to_origin(Surface_mesh &mesh) {
    BBox bbox = CGAL::bounding_box(mesh.points().begin(), mesh.points().end());

    FT x = (bbox.xmin() + bbox.xmax()) / 2;
    FT y = (bbox.ymin() + bbox.ymax()) / 2;
    FT z = (bbox.zmin() + bbox.zmax()) / 2;
    Vector center = Vector(x, y, z);

    for (Vertex_index v : vertices(mesh)) {
        Mesh_point p = mesh.point(v);
        mesh.point(v) = Mesh_point(p.x() - center.x(), p.y() - center.y(),
                                   p.z() - center.z());
    }
    return center;
}

void move_mesh(Surface_mesh &mesh, Vector center) {
    for (Vertex_index v : vertices(mesh)) {
        Mesh_point p = mesh.point(v);
        mesh.point(v) = Mesh_point(p.x() - center.x(), p.y() - center.y(),
                                   p.z() - center.z());
    }
}

/***********************Geometry***************************/

Plane* normalize_plane(const Plane* raw_plane)
{
    Eigen::Vector3f normal (CGAL::to_double(raw_plane->a()),
                            CGAL::to_double(raw_plane->b()),
                            CGAL::to_double(raw_plane->c()));
    float norm = normal.norm();
    normal /= norm;
    float d = CGAL::to_double(raw_plane->d()) / norm;
    return new Plane(normal(0), normal(1), normal(2), d);
}

void calculate_areas_faces_3d(Surface_mesh &mesh) {

    bool created;
    Surface_mesh::Property_map<Face_index, float> face_areas;
    boost::tie(face_areas, created) =
        mesh.add_property_map<Face_index, float>("f:area", 0.0);
    assert(created);

    Surface_mesh::Property_map<Face_index, Vector> face_normal;
    boost::tie(face_normal, created) =
        mesh.add_property_map<Face_index, Vector>("f:normal_", Vector(0, 0, 0));
    assert(created);
    PMP::compute_face_normals(mesh, face_normal);

    for (Face_index f : faces(mesh)) {
        Eigen::Vector3d unit_normal(CGAL::to_double(face_normal[f].x()),
                                    CGAL::to_double(face_normal[f].y()),
                                    CGAL::to_double(face_normal[f].z()));
        unit_normal.normalize();
        Eigen::Vector3d total(0, 0, 0);
        Halfedge_index e = mesh.halfedge(f);
        Halfedge_index begin_e = e;
        do {
            Vertex_index v = mesh.source(e);
            Vertex_index v_next = mesh.target(e);
            Mesh_point p = mesh.point(v);
            Mesh_point p_next = mesh.point(v_next);
            Eigen::Vector3d v1(CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                               CGAL::to_double(p.z()));
            Eigen::Vector3d v2(CGAL::to_double(p_next.x()),
                               CGAL::to_double(p_next.y()),
                               CGAL::to_double(p_next.z()));
            total += v1.cross(v2);

            e = mesh.next(e);
        } while (e != begin_e);

        face_areas[f] = std::abs(total.dot(unit_normal) / 2);
    }
}

float calculate_area_cell_2(const Dart_descriptor &face_dart, LCC_3 &lcc,
                            const Plane *plane) {
    assert(plane != nullptr);
    Vector normal_ = plane->orthogonal_vector();
    Eigen::Vector3f normal(CGAL::to_double(normal_.x()),
                           CGAL::to_double(normal_.y()),
                           CGAL::to_double(normal_.z()));
    normal.normalize();
    Eigen::Vector3f total(0, 0, 0);

    Dart_descriptor v_dart = face_dart;
    do {
        Dart_descriptor next_v_dart = lcc.beta<1>(v_dart);
        Eigen::Vector3f v1(CGAL::to_double(lcc.point(v_dart).x()),
                           CGAL::to_double(lcc.point(v_dart).y()),
                           CGAL::to_double(lcc.point(v_dart).z()));
        Eigen::Vector3f v2(CGAL::to_double(lcc.point(next_v_dart).x()),
                           CGAL::to_double(lcc.point(next_v_dart).y()),
                           CGAL::to_double(lcc.point(next_v_dart).z()));
        total += v1.cross(v2);
        v_dart = next_v_dart;
    } while (v_dart != face_dart);

    return std::abs(total.dot(normal)) / 2;
}

float calculate_area_cell_2(
    const Dart_descriptor &face_dart, LCC_3 &lcc) {
    // calculate the normal
    const Plane *plane = lcc.info<2>(face_dart).second;
    float area = calculate_area_cell_2(face_dart, lcc, plane);
    return area;
}

float dist_pt_to_plane(const Point &pt, const Plane &plane) {
    float dist =
        std::abs(CGAL::to_double(plane.a() * pt.x() + plane.b() * pt.y() +
                                 plane.c() * pt.z() + plane.d()) /
                 std::sqrt(CGAL::to_double(plane.a() * plane.a() +
                                           plane.b() * plane.b() +
                                           plane.c() * plane.c())));
    return dist;
}

bool almost_same_plane(
    const BBox &bbox,
    const Plane *a, const Plane *b,
    float angle_thresh = 10,
    float dist_thresh = 0.5)
{
     Vector an = a->orthogonal_vector(), bn = b->orthogonal_vector();
    float angle = std::acos(CGAL::to_double(an * bn) /
                            std::sqrt(CGAL::to_double(an.squared_length() *
                                                      bn.squared_length()))) *
                  180 / M_PI;
    if (angle > angle_thresh && angle < 180 - angle_thresh) // not parallel
        return false;
    // if (!CGAL::parallel(*a, *b)) return false;

    // TODO: ??? close distance to origin along the orthogonal vector
    // FT ka = -a->d() / an.squared_length(), kb = -b->d() / bn.squared_length();
    // Point pa(0, 0, 0), pb(0, 0, 0);
    // pa += an * ka;
    // pb += bn * kb;
    // if (CGAL::squared_distance(pa, pb) > dist_thresh * dist_thresh)
    //     return false;

    // return true;

    FT xmin = bbox.xmin(), xmax = bbox.xmax(),
       ymin = bbox.ymin(), ymax = bbox.ymax(),
       zmin = bbox.zmin(), zmax = bbox.zmax();
    Point p1(xmin, ymin, zmin), p2(xmax, ymin, zmin), p3(xmax, ymax, zmin),
          p4(xmin, ymax, zmin), p5(xmin, ymin, zmax), p6(xmax, ymin, zmax),
          p7(xmax, ymax, zmax), p8(xmin, ymax, zmax);
    
    std::vector<Point> proj_a {
            a->projection(p1), a->projection(p2), a->projection(p3), a->projection(p4),
            a->projection(p5), a->projection(p6), a->projection(p7), a->projection(p8)};
    std::vector<float> dists {
        dist_pt_to_plane(proj_a[0], *b), dist_pt_to_plane(proj_a[1], *b), dist_pt_to_plane(proj_a[2], *b), dist_pt_to_plane(proj_a[3], *b),
        dist_pt_to_plane(proj_a[4], *b), dist_pt_to_plane(proj_a[5], *b), dist_pt_to_plane(proj_a[6], *b), dist_pt_to_plane(proj_a[7], *b)};
    
    float max_dist = *std::max_element(dists.begin(), dists.end());
    // float ave_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    if (max_dist > dist_thresh)
        return false;
    return true;    
}

float dist_pt_to_plane(
    const Point &pt, const Plane *plane)
{
    float dist = std::abs(CGAL::to_double(plane->a() * pt.x() + plane->b() * pt.y() +
                                         plane->c() * pt.z() + plane->d()) /
                          std::sqrt(CGAL::to_double(plane->a() * plane->a() +
                                                    plane->b() * plane->b() +
                                                    plane->c() * plane->c())));
    return dist;
}

bool almost_same_plane(
    const Surface_mesh &mesh,
    const Plane *a, const std::vector<Face_index> &faces_a,
    const Plane *b, const std::vector<Face_index> &faces_b,
    float angle_thresh = 1,
    float dist_thresh = 0.1) 
{
    
    Eigen::Vector3d an(CGAL::to_double(a->a()), CGAL::to_double(a->b()), CGAL::to_double(a->c()));
    Eigen::Vector3d bn(CGAL::to_double(b->a()), CGAL::to_double(b->b()), CGAL::to_double(b->c()));
    an.normalize();
    bn.normalize();
    float radian = std::acos(an.dot(bn));
    if (radian > angle_thresh/180 * M_PI && radian < M_PI - angle_thresh/180*M_PI)
    {
        // std::cout << "radian: " << radian << ", angle_thresh: " << angle_thresh/180*M_PI << '\n';
        return false;
    }
        
    
    float max_dist_a_to_b = 0, max_dist_b_to_a = 0;
    for (Face_index f : faces_a) {
        for (Halfedge_index e : halfedges_around_face(halfedge(f, mesh), mesh)) {
            Vertex_index v = target(e, mesh);
            Point p = mesh.point(v);
            float dist = dist_pt_to_plane(p, b);
            max_dist_a_to_b = std::max(max_dist_a_to_b, dist);
        }
    }

    for (Face_index f : faces_b) {
        for (Halfedge_index e : halfedges_around_face(halfedge(f, mesh), mesh)) {
            Vertex_index v = target(e, mesh);
            Point p = mesh.point(v);
            float dist = dist_pt_to_plane(p, a);
            max_dist_b_to_a = std::max(max_dist_b_to_a, dist);
        }
    }

    if (max_dist_a_to_b > dist_thresh
    || max_dist_b_to_a > dist_thresh)
        return false;
    else
        return true;
}

std::pair<bool, const Plane *>
is_duplicate_plane(
    const Surface_mesh &mesh,
    const Plane *to_check_plane, const std::vector<Face_index> &to_check_faces,
    const std::unordered_map<const Plane*, std::vector<Face_index>> &plane_regions,
    float angle_thresh = 10, 
    float dist_thresh = 0.5) {

    for (const auto &[plane, faces] : plane_regions) {
        if (to_check_plane == plane)
            return std::make_pair(true, plane);

        if (almost_same_plane(mesh, to_check_plane, to_check_faces, plane, faces, angle_thresh, dist_thresh))
            return std::make_pair(true, plane);
    }

    return std::make_pair(false, nullptr);
}

std::pair<bool, const Plane *>
is_duplicate_plane(
    const BBox &bbox,
    const Plane *to_check_plane,
    const std::vector<const Plane *> &planes,
    float angle_thresh = 1, float dist_thresh = 0.1) {
    for (const auto &plane : planes) {
        if (to_check_plane == plane)
            return std::make_pair(true, plane);

        if (almost_same_plane(bbox, to_check_plane, plane, angle_thresh, dist_thresh))
            return std::make_pair(true, plane);
    }
    return std::make_pair(false, nullptr);
}

std::pair<bool, const Plane *>
is_duplicate_plane(
    const BBox &bbox,
    const Plane *to_check_plane,
    const std::set<const Plane *> &planes,
    float angle_thresh = 1, float dist_thresh = 0.1) {
    for (const auto &plane : planes) {
        if (to_check_plane == plane)
            return std::make_pair(true, plane);

        if (almost_same_plane(bbox, to_check_plane, plane, angle_thresh, dist_thresh))
            return std::make_pair(true, plane);
    }
    return std::make_pair(false, nullptr);
}

float radian_between_two_faces(const Dart_descriptor &this_face,
                               const Dart_descriptor &next_face, LCC_3 &lcc,
                               const Plane *this_plane,
                               const Plane *next_plane) {
    assert(this_plane != nullptr && next_plane != nullptr);
    if (this_plane == next_plane)
        return M_PI;

    Vector this_normal = this_plane->orthogonal_vector(),
           next_normal = next_plane->orthogonal_vector();
    float dot = CGAL::to_double(this_normal * next_normal);
    float norm_this = std::sqrt(CGAL::to_double(this_normal.squared_length())),
          norm_next = std::sqrt(CGAL::to_double(next_normal.squared_length()));
    float radian = std::acos(dot / (norm_this * norm_next)); // from 0 to pi
    // std::cout << "cos radian: " << radian << '\n';

    Dart_descriptor vtx_on_this_face = lcc.previous(this_face),
                    vtx_on_next_face = lcc.previous(next_face);
    Point pt_on_this_face = lcc.point(vtx_on_this_face),
          pt_on_next_face = lcc.point(vtx_on_next_face);
    CGAL::Oriented_side this_side, next_side;
    this_side = this_plane->oriented_side(pt_on_next_face);
    if (this_side == CGAL::ON_ORIENTED_BOUNDARY) {
        for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator
                 it(lcc.one_dart_per_incident_cell<0, 3>(next_face).begin()),
             itend(lcc.one_dart_per_incident_cell<0, 3>(next_face).end());
             it != itend; ++it) {
            this_side = this_plane->oriented_side(lcc.point(it));
            if (this_side != CGAL::ON_ORIENTED_BOUNDARY)
                break;
        }
    }
    // while(this_side == CGAL::ON_ORIENTED_BOUNDARY && vtx_on_next_face !=
    // next_face)
    // {
    //     vtx_on_next_face = lcc.previous(vtx_on_next_face);
    //     pt_on_next_face = lcc.point(vtx_on_next_face);
    //     this_side = this_plane->oriented_side(pt_on_next_face);
    // }
    if (this_side == CGAL::ON_ORIENTED_BOUNDARY) {
        //     std::cout << "this_side: " << this_side << '\n';
        //     std::cout << "this_plane: " << this_plane << " "<<
        //     this_plane->a() << " " << this_plane->b() << " " <<
        //     this_plane->c() << " " << this_plane->d() << '\n'; std::cout <<
        //     "next_plane: " << next_plane << " "<< next_plane->a() << " " <<
        //     next_plane->b() << " " << next_plane->c() << " " <<
        //     next_plane->d() << '\n';
        return M_PI;
    }
    next_side = next_plane->oriented_side(pt_on_this_face);
    if (next_side == CGAL::ON_ORIENTED_BOUNDARY) {
        for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator
                 it(lcc.one_dart_per_incident_cell<0, 3>(this_face).begin()),
             itend(lcc.one_dart_per_incident_cell<0, 3>(this_face).end());
             it != itend; ++it) {
            next_side = next_plane->oriented_side(lcc.point(it));
            if (next_side != CGAL::ON_ORIENTED_BOUNDARY)
                break;
        }
    }
    // while(next_side == CGAL::ON_ORIENTED_BOUNDARY && vtx_on_this_face !=
    // this_face)
    // {
    //     vtx_on_this_face = lcc.previous(vtx_on_this_face);
    //     pt_on_this_face = lcc.point(vtx_on_this_face);
    //     next_side = next_plane->oriented_side(pt_on_this_face);
    // }
    if (next_side == CGAL::ON_ORIENTED_BOUNDARY) {
        // std::cout << "next_side: " << next_side << '\n';
        // std::cout << "this_plane: " << this_plane << " "<< this_plane->a() <<
        // " " << this_plane->b() << " " << this_plane->c() << " " <<
        // this_plane->d() << '\n'; std::cout << "next_plane: " << next_plane <<
        // " "<< next_plane->a() << " " << next_plane->b() << " " <<
        // next_plane->c() << " " << next_plane->d() << '\n';
        return M_PI;
    }
    // std::cout << " side: " << this_side << " " << next_side << '\n';
    if (this_side == next_side)
        radian = M_PI - radian;
    return radian;
}

int count_in_cell(
    std::unordered_map<int, std::pair<bool, Dart_descriptor>> is_in_cell) {
    int num_in_cell = 0;
    for (const auto &[cell, is_in] : is_in_cell)
        if (is_in.first)
            num_in_cell++;
    return num_in_cell;
}

Point2 center_of_convex_hull(const std::vector<Point2> &ch_pts) {
    // https://graphics.stanford.edu/courses/cs368-04-spring/manuals/CGAL_Tutorial.pdf,
    // Pages 18-19 Origin(0, 0)
    assert(ch_pts.size() > 2);
    int n = ch_pts.size();
    Vector2 center(0, 0);
    FT denominator = 0;
    for (int i = 0; i < n; ++i) {
        int cur = i;
        int next = (i + 1) % n;
        FT a = ch_pts[cur].x() * ch_pts[next].y() -
               ch_pts[next].x() * ch_pts[cur].y();
        center = center + a * (Vector2(ch_pts[cur].x(), ch_pts[cur].y()) +
                               Vector2(ch_pts[next].x(), ch_pts[cur].y()));
        denominator += a;
    }
    center = center / (3 * denominator);
    // std::cout << "center of convex hull 2: " << center << std::endl;
    return Point2(center.x(), center.y());
}

void extend_convex_hull(std::vector<Point2> &ch_pts, const float &ext_dist) {
    Point2 center = center_of_convex_hull(ch_pts);
    for (auto &p : ch_pts) {
        // std::cout << "p: " << p;
        Vector2 v = p - center;
        p = p + 0.1 * v;
        // std::cout << " " << p << "\n";
    }
}

void compute_extended_convex_hull_of_planes(
    const std::unordered_map<const Plane *, std::vector<Face_index>>
        plane_regions,
    std::unordered_map<const Plane *, const std::vector<Point>>
        &plane_extended_convex_hulls,
    Surface_mesh &input_mesh, const float &ext_dist = 0.1) {
    for (const auto &[plane, faces] : plane_regions) {
        // std::cout << "plane: " << plane << '\n';
        std::vector<Point2> pts;
        for (const auto &f : faces) {
            Halfedge_index e = input_mesh.halfedge(f);
            Halfedge_index begin_e = e;
            do {
                Mesh_point p = input_mesh.point(input_mesh.source(e));
                // Point p3 = Point(p.x(), p.y(), p.z());
                Point2 p2 = plane->to_2d(p);
                // Point2 p2 = project_pt3d_on_plane2d(p3, *plane);
                pts.push_back(p2);
                e = input_mesh.next(e);
            } while (e != begin_e);
        }
        // get the 2d convex hull
        std::vector<Point2> ch_pts;
        CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(ch_pts));
        // extend the extreme points of the convex hull by a distance
        extend_convex_hull(ch_pts, ext_dist);
        // back project the 2d points to 3d points by plane->to_3d(p2)
        std::vector<Point> ech_pts;
        for (const auto &p2 : ch_pts) {
            Point p3 = plane->to_3d(p2);
            // Point p3 = backproject(*plane, p2, 0);
            ech_pts.push_back(p3);
        }

        plane_extended_convex_hulls.insert({plane, ech_pts});
    }
}

/***********************VIS***************************/

// void draw_sorting(
//     const std::vector<std::pair<Dart_descriptor, std::pair<float, Eva_cut>>>
//     top_planes_of_cells, const std::unordered_map<Dart_descriptor,
//     std::vector<std::pair<const Plane*, float>>> to_cut_planes, LCC_3 &lcc,
//     Surface_mesh &input_mesh,
//     const std::unordered_map<const Plane*, std::vector<Face_index>>
//     &plane_regions)
// {
//     int num_cell = top_planes_of_cells.size();
//     bool found;
//     Surface_mesh::Property_map<Face_index, CGAL::IO::Color> f_colors;
//     boost::tie(f_colors, found) =
//         input_mesh.property_map<Face_index, CGAL::IO::Color>("f:color");

//     std::cout << "enter the cell id to draw(0-" << num_cell - 1 << "): ";
//     int cid = 0;

//     // while(std::cin >> cid)
//     // {
//     //     if (cid >= 0 && cid < num_cell)
//     //     {
//             // colorize the input mesh as gray
//             for (auto f : input_mesh.faces())
//                 f_colors[f] = CGAL::IO::Color(128, 128, 128);

//             Dart_descriptor cell = top_planes_of_cells[cid].first;
//             std::vector<std::pair<const Plane*, float>> planes =
//             to_cut_planes.at(cell);
//             // colorize the regions of the top planes as red; the regions of
//             the other planes as green const Plane* top_plane =
//             planes[0].first; std::vector<Face_index> top_region =
//             plane_regions.at(top_plane); for (const auto f : top_region)
//                 f_colors[f] = CGAL::IO::Color(255, 0, 0);

//             for (std::vector<std::pair<const Plane*, float>>::const_iterator
//                     p = planes.begin() + 1; p != planes.end(); ++p)
//             {
//                 std::vector<Face_index> region =
//                 plane_regions.at((*p).first); for (const auto f : region)
//                     f_colors[f] = CGAL::IO::Color(0, 255, 0);
//             }

//             Surface_mesh to_draw_mesh = input_mesh;
//             Surface_mesh::Property_map<Face_index, CGAL::IO::Color>
//             df_colors; boost::tie(df_colors, found) =
//                 to_draw_mesh.property_map<Face_index,
//                 CGAL::IO::Color>("f:color");
//             // add vertices of the cell in the to_draw mesh
//             // add faces of the cell in the to_draw mesh
//             LCC_3::One_dart_per_incident_cell_range<2, 3> faces =
//             lcc.one_dart_per_incident_cell<2, 3>(cell); for
//             (LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator
//                 it=faces.begin(); it!=faces.end(); ++it) // the face of lcc,
//                 may not be a face in surface mesh due to precision issue
//             {
//                 Dart_descriptor d = it;
//                 std::vector<Vertex_index> vidxs;
//                 do
//                 {
//                     Point pt (lcc.point(d).x(), lcc.point(d).y(),
//                     lcc.point(d).z()); Vertex_index vid =
//                     to_draw_mesh.add_vertex(pt); vidxs.push_back(vid); d =
//                     lcc.beta<1>(d);
//                 }while(d != it);
//                 Face_index f = to_draw_mesh.add_face(vidxs);
//                 df_colors[f] = CGAL::IO::Color(0, 0, 255);
//             }
//             PMP::triangulate_faces(to_draw_mesh);
//             // PMP::remove_almost_degenerate_faces(to_draw_mesh);
//             CGAL::draw(to_draw_mesh); // bug: degenerate faces in the to draw
//             mesh std::cout << "enter the cell id to draw(0-" << num_cell - 1
//             << "): ";
//     //     }
//     //     else
//     //     {
//     //         break;
//     //     }
//     // }

// }

float length_of_cell1(
    const LCC_3 &lcc,
    Dart_const_descriptor &edge)
{   
    FT l = CGAL::squared_distance(
            lcc.point(edge), lcc.point(lcc.beta<1>(edge)));;
    float length = std::sqrt(CGAL::to_double(l));

    return length;
}

float area_of_cell2(
    const LCC_3 &lcc,
    Dart_const_descriptor face)
{
    FT area = 0;
    const Plane *plane = lcc.info<2>(face).second;
    std::vector<Point2> pts;
    Dart_const_descriptor d = face;
    do
    {
        pts.push_back(plane->to_2d(lcc.point(d)));
        d = lcc.beta<1>(d);
    } while (d != face);

    for (int i = 0; i < pts.size(); i++)
    {
        FT x1 = pts[i].x(), 
           y1 = pts[i].y();
        FT x2 = pts[(i+1)%pts.size()].x(), 
           y2 = pts[(i+1)%pts.size()].y();
        area += (x1*y2 - x2*y1) / 2;
    }
    return std::abs(CGAL::to_double(area));
}

float volume_of_cell3(
    const LCC_3 &lcc,
    Dart_const_descriptor cell)
{
    FT volume = 0;

    for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator 
        f = lcc.one_dart_per_incident_cell<2, 3>(cell).begin();
        f != lcc.one_dart_per_incident_cell<2, 3>(cell).end(); ++f)
    {
        Dart_const_descriptor d = f;
        std::vector<Point> pts;
        do
        {
            pts.push_back(lcc.point(d));
            d = lcc.beta<1>(d);
        } while (d != f);

        // normal of triangles pointing inward of the cell
        FT x1 = pts[0].x(), y1 = pts[0].y(), z1 = pts[0].z();
        for (int i = 1; i < pts.size()-1; i++)
        {
            FT x2 = pts[i].x(), y2 = pts[i].y(), z2 = pts[i].z();
            FT x3 = pts[i+1].x(), y3 = pts[i+1].y(), z3 = pts[i+1].z();
            FT v = 1.0/6.0 * 
                    (-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
            volume += v;
        }
    }
    return std::abs(CGAL::to_double(volume));
}   

float volume_of_cell3(
    const LCC_3 &lcc,
    const Volume_idx &vid)
{
   assert (vid >= 0 && vid < lcc.number_of_attributes<3>());
   
   for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (lcc.info<3>(it).first == vid)
            return volume_of_cell3(lcc, it);
    } 
}

float volume_of_almost_closed_surface_mesh(
    const Surface_mesh &mesh)
{
    FT volume = 0;
    for (Face_index f : faces(mesh))
    {
        std::vector<Point> pts;
        Halfedge_index e = mesh.halfedge(f);
        Halfedge_index begin_e = e;
        do
        {
            pts.push_back(mesh.point(mesh.source(e)));
            e = mesh.next(e);
        } while (e != begin_e);

        // normal of triangles pointing inward of the cell
        FT x1 = pts[0].x(), y1 = pts[0].y(), z1 = pts[0].z();
        for (int i = 1; i < pts.size()-1; i++)
        {
            FT x2 = pts[i].x(), y2 = pts[i].y(), z2 = pts[i].z();
            FT x3 = pts[i+1].x(), y3 = pts[i+1].y(), z3 = pts[i+1].z();
            FT v = 1.0/6.0 * 
                    (-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
            volume += v;
        }
    }

    return std::abs(CGAL::to_double(volume));
}


float area_of_cell3(
    const LCC_3 &lcc,
    Dart_const_descriptor &cell)
{
    float area = 0;
    for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator
        f = lcc.one_dart_per_incident_cell<2, 3>(cell).begin();
        f != lcc.one_dart_per_incident_cell<2, 3>(cell).end(); ++f)
    {
        area += area_of_cell2(lcc, f);
    }

    return area;
}

float area_of_cell3(
    const LCC_3 &lcc,
    const Volume_idx &vid)
{
    assert (vid >= 0 && vid < lcc.number_of_attributes<3>());
    
    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (lcc.info<3>(it).first == vid)
            return area_of_cell3(lcc, it);
    }
}

float convex_hull_volume_of_cell3(
    const LCC_3 &lcc,
    Dart_const_descriptor &cell,
    Surface_mesh &ch_mesh)
{
    std::vector<Point> pts;
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::const_iterator v =
             lcc.one_dart_per_incident_cell<0, 3>(cell).begin();
         v != lcc.one_dart_per_incident_cell<0, 3>(cell).end(); ++v)
        pts.push_back(lcc.point(v));

    CGAL::convex_hull_3(pts.begin(), pts.end(), ch_mesh);

    if (!CGAL::is_closed(ch_mesh))
        return 0;

    return CGAL::to_double(PMP::volume(ch_mesh));
}

float convex_hull_volume_of_cell3(
    const LCC_3 &lcc,
    const Volume_idx &vid,
    Surface_mesh &ch_mesh)
{
    assert (vid >= 0 && vid < lcc.number_of_attributes<3>());
    
    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (lcc.info<3>(it).first == vid)
            return convex_hull_volume_of_cell3(lcc, it, ch_mesh);
    }
}

void clean_lcc_edges(
    LCC_3 &lcc)
{
    std::vector<Dart_descriptor> edges_to_remove;
    for (LCC_3::One_dart_per_cell_range<1>::iterator
        it=lcc.one_dart_per_cell<1>().begin(),
        itend=lcc.one_dart_per_cell<1>().end(); it!=itend; ++it)
    {
        std::set<const Plane*> planes;
        for (LCC_3::One_dart_per_incident_cell_range<2, 1>::iterator
            f=lcc.one_dart_per_incident_cell<2, 1>(it).begin(),
            fend=lcc.one_dart_per_incident_cell<2, 1>(it).end(); f!=fend; ++f)
            planes.insert(lcc.info<2>(f).second);
        if (planes.size() == 1)
            edges_to_remove.push_back(it);
    }

    for (auto edge : edges_to_remove)
        lcc.remove_cell<1>(edge);
    
    // re-index the faces
    Vertex_idx vtx_id = 0;
    for (LCC_3::One_dart_per_cell_range<0>::iterator
        it=lcc.one_dart_per_cell<0>().begin(),
        itend=lcc.one_dart_per_cell<0>().end(); it!=itend; ++it)
        lcc.info<0>(it) = vtx_id++;

    Face_idx fid = 0;
    for (LCC_3::One_dart_per_cell_range<2>::iterator
        it=lcc.one_dart_per_cell<2>().begin(),
        itend=lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
        lcc.info<2>(it).first = fid++;

    Volume_idx vid = 0;
    for (LCC_3::One_dart_per_cell_range<3>::iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        lcc.info<3>(it).first = vid++;
        lcc.info<3>(it).second = true;
    }
}

/***********************IO***************************/

void editFirstLine(const std::string &fileName, std::size_t n)
{
    if (n == 0)
        return;
    std::vector<std::string> lines;
    std::string _line;
    std::ifstream inputFile(fileName);
    while (std::getline(inputFile, _line))
    {
        lines.push_back(_line);
    }
    inputFile.close();
    if(lines.size() < 2) return;

    std::vector<std::string> results;
    boost::split(results, lines[1], boost::is_any_of(" "));
    if(results.size() != 3) return;
    results[1] = std::to_string(n);
    lines[1] = boost::algorithm::join(results, " ");

    std::ofstream outputFile(fileName, std::ios::trunc);
    for (const auto &l : lines)
    {
        outputFile << l << "\n";
    }
    outputFile.close();
}

bool write_lcc_as_off(
    const LCC_3 &i_lcc, 
    std::string off_file, 
    bool in_only, 
    bool write_sta = false,
    std::string sta_file = "")
{
    LCC_3 lcc = i_lcc;
    // modified from
    // https://github.com/CGAL/cgal/blob/ed3503d2381842f837a8118a749ab380f889485b/Linear_cell_complex/include/CGAL/Linear_cell_complex_constructors.h#L299
    std::ofstream out(off_file.c_str());

    if (!lcc.are_all_faces_closed()) {
        std::cerr << "Impossible to write in off a map having open faces."
                  << std::endl;
        return false;
    }

    Vertex_idx vtx_id = 0;
    for (LCC_3::One_dart_per_cell_range<0>::iterator
        it=lcc.one_dart_per_cell<0>().begin(),
        itend=lcc.one_dart_per_cell<0>().end(); it!=itend; ++it)
        lcc.info<0>(it) = vtx_id++;

    CGAL::File_header_OFF header(false);
    // header.set_binary(IO::is_binary(out));
    header.set_no_comments(!CGAL::IO::is_pretty(out));
    CGAL::File_writer_OFF writer(header);
    // writer.header().set_polyhedral_surface(true);
    writer.header().set_halfedges(lcc.number_of_darts());
    writer.header().set_verbose(true);

    // Print header.
    writer.write_header(out, lcc.number_of_vertex_attributes(),
                        lcc.number_of_halfedges(),
                        lcc.template one_dart_per_cell<2>().size() * 2);

    typedef typename LCC_3::Vertex_attribute_range::iterator VCI;
    VCI vit, vend = lcc.vertex_attributes().end();
    for (vit = lcc.vertex_attributes().begin(); vit != vend; ++vit) {
        writer.write_vertex(::CGAL::to_double(vit->point().x()),
                            ::CGAL::to_double(vit->point().y()),
                            ::CGAL::to_double(vit->point().z()));
    }

    typedef CGAL::Inverse_index<VCI> Index;
    Index index(lcc.vertex_attributes().begin(), lcc.vertex_attributes().end());
    writer.write_facet_header();
    
    int num_cell3 = lcc.number_of_attributes<3>();
    std::vector<CGAL::IO::Color> colors;
    create_colormap(num_cell3, colors);
    std::vector<bool> sta_cell3(num_cell3, true);

    std::size_t cnt = 0;
    typename LCC_3::size_type m = lcc.get_new_mark();
    for (typename LCC_3::Dart_range::iterator itall = lcc.darts().begin(),
                                              itallend = lcc.darts().end();
         itall != itallend; ++itall) {
        if (!lcc.is_marked(itall, m)) {

            if (in_only && !lcc.info<3>(itall).second)
            {
                sta_cell3[lcc.info<3>(itall).first] = false;
                continue;
            }

            std::size_t n = 0;
            typename LCC_3::Dart_handle cur = itall;
            do {
                ++n;
                assert(lcc.is_next_exist(cur));
                cur = lcc.next(cur);
            } while (cur != itall);

            CGAL_assertion(n >= 3);
            writer.write_facet_begin(n);

            // Second we write the indices of vertices.
            do {
                writer.write_facet_vertex_index(
                    index[VCI(lcc.vertex_attribute(cur))]);

                lcc.mark(cur, m);
                lcc.mark(lcc.other_orientation(cur), m);                    // for GMap only, for CMap  marks the same dart twice
                assert(lcc.is_next_exist(cur)); 
                cur = lcc.next(cur);
            } while (cur != itall);

            writer.write_face_color(
                colors[lcc.info<3>(itall).first].red(),
                colors[lcc.info<3>(itall).first].green(),
                colors[lcc.info<3>(itall).first].blue());

            writer.write_facet_end();
            cnt++;
        }
    }

    // for (typename LCC_3::One_dart_per_cell_range<3>::iterator 
    //     it=lcc.one_dart_per_cell<3>().begin(), itend=lcc.one_dart_per_cell<3>().end();
    //     it != itend; ++it)
    // {
    //     if (in_only && !lcc.info<3>(it).second)
    //     {
    //         sta_cell3[lcc.info<3>(it).first] = false;
    //         continue;
    //     }

    //     for (typename LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator 
    //         f=lcc.one_dart_per_incident_cell<2, 3>(it).begin(),
    //         fend=lcc.one_dart_per_incident_cell<2, 3>(it).end(); f!=fend; ++f)
    //     {
    //         size_t n = 0;
    //         typename LCC_3::Dart_handle cur = f;
    //         do
    //         {
    //             ++n;
    //             assert(lcc.is_next_exist(cur));
    //             cur = lcc.next(cur);
    //         } while (cur != f);

    //         CGAL_assertion(n >= 3);
    //         writer.write_facet_begin(n);
            
    //         do
    //         {
    //             writer.write_facet_vertex_index(
    //                 index[VCI(lcc.vertex_attribute(cur))]);
    //             assert(lcc.is_next_exist(cur)); 
    //             cur = lcc.next(cur);
    //         } while (cur != f);

    //         writer.write_face_color(
    //             colors[lcc.info<3>(it).first].red(),
    //             colors[lcc.info<3>(it).first].green(),
    //             colors[lcc.info<3>(it).first].blue());
            
    //         writer.write_facet_end();
    //     }
    // }

    writer.write_footer();
    // lcc.free_mark(m);
    out.close();
    editFirstLine(off_file, cnt);
    std::cout << "save " << off_file << " successfully\n";
    
    if (write_sta)
    {
        std::ofstream ofs(sta_file);
        for (int i = 0; i < num_cell3; i++)
        {
            if (sta_cell3[i])
            {
                // std::cout << "calculating volume: " << i << "\n";
                float volume = volume_of_cell3(lcc, i);
                // std::cout << "\t vol = " << volume
                //         << " r=" << static_cast<int>(colors[i].red()) 
                //         << " g=" << static_cast<int>(colors[i].green())
                //         << " b=" << static_cast<int>(colors[i].blue()) << "\n";
                ofs << "\t vol = " << volume
                    << " r=" << static_cast<int>(colors[i].red()) 
                    << " g=" << static_cast<int>(colors[i].green())
                    << " b=" << static_cast<int>(colors[i].blue()) << "\n";
            }
        }
        ofs.close();
        std::cout << "save " << sta_file << " successfully\n";
    }
    return true;
}


struct Evaluate_weights {
    float volume = 1.0;
    float ratio = 1.0;
    float nxy = 1.0;
    float nz = 1.0;
    float adjacency = 1.0;
};

enum DECOM {
    CONVEXITY,
    VOLUME_LARGE_STD,
    VOLUME_SMALL_STD,
    RATIO_LARGE_STD,
    RATIO_SMALL_STD,
    ADJACENCY_LARGE_TOP,
    ADJACENCY_SMALL_TOP,
    NXY,
    VERTICAL,
    HORIZONTAL,
    CENTRAL,
    PERIPHERAL,
    // CENTRAL_VERTICAL,
    // CENTRAL_HORIZONTAL
};

float data_mean(std::vector<float> &data) {
    float sum = 0.0, mean;
    int i;

    for (i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    mean = sum / data.size();
    return mean;
}

float data_std(std::vector<float> &data) {
    float mean, standardDeviation = 0.0;

    mean = data_mean(data);
    int i;
    for (i = 0; i < data.size(); ++i) {
        standardDeviation += std::pow(data[i] - mean, 2);
    }

    return std::sqrt(standardDeviation / data.size());
}

float get_pitch(const Plane *p) {
    Vector _n = p->orthogonal_vector();
    Eigen::Vector3f n(CGAL::to_double(_n.x()), CGAL::to_double(_n.y()),
                      CGAL::to_double(_n.z()));
    n.normalize();
    float radian =
        std::acos(std::sqrt(1 - n(2) * n(2))); // [0 vertical, pi/2 horizontal]
    return radian;
}

void dfs(int node, const std::vector<std::vector<int>> &adjList,
         std::vector<bool> &visited, std::set<int> &component) {
    visited[node] = true;
    component.insert(node);

    for (int neighbor : adjList[node]) {
        if (!visited[neighbor]) {
            dfs(neighbor, adjList, visited, component);
        }
    }
}

void connected_component(const std::vector<std::pair<int, int>> &edges,
                         std::vector<std::set<int>> &components) {
    int numNodes = 0;
    std::unordered_map<int, int> nodeMap;
    std::unordered_map<int, int> invNodeMap;
    std::vector<std::vector<int>> adjList;

    for (const auto &edge : edges) {
        int u = edge.first;
        int v = edge.second;

        if (nodeMap.find(u) == nodeMap.end()) {
            invNodeMap[numNodes] = u;
            nodeMap[u] = numNodes++;
            adjList.push_back({});
        }
        if (nodeMap.find(v) == nodeMap.end()) {
            invNodeMap[numNodes] = v;
            nodeMap[v] = numNodes++;
            adjList.push_back({});
        }

        int uIndex = nodeMap[u];
        int vIndex = nodeMap[v];

        adjList[uIndex].push_back(vIndex);
        adjList[vIndex].push_back(uIndex);
    }

    std::vector<bool> visited(numNodes, false);
    for (int i = 0; i < numNodes; i++) {
        if (!visited[i]) {
            std::set<int> temp_comp;
            dfs(i, adjList, visited, temp_comp);
            // std::cout << "found comps: ";
            std::set<int> component;
            for (auto &c : temp_comp) {
                component.insert(invNodeMap[c]);
                // std::cout << invNodeMap[c] << " ";
            }
            components.push_back(component);
        }
    }
}

Eigen::Vector2f get_xy_normal_of_vertical_plane(const Plane *p) {
    Vector _n = p->orthogonal_vector();
    Eigen::Vector2f n_xy(CGAL::to_double(_n.x()), CGAL::to_double(_n.y()));
    n_xy.normalize();
    return n_xy;
}

bool are_vertical_planes_orthogonal(Eigen::Vector2f &n1, Eigen::Vector2f &n2,
                                    const float &thresh = 0.2) {
    n1.normalize();
    n2.normalize();
    double dot = n1.dot(n2);
    if (std::abs(dot) < thresh)
        return true;
    else
        return false;
}

bool are_vertical_planes_orthogonal(
    const Plane *p1, const Plane *p2, const float &thresh = 0.2) {
    Eigen::Vector2f n1 = get_xy_normal_of_vertical_plane(p1);
    Eigen::Vector2f n2 = get_xy_normal_of_vertical_plane(p2);
    return are_vertical_planes_orthogonal(n1, n2, thresh);
}

bool is_line_horizontal(Eigen::Vector3f &n, float thresh = 0.2) {
    n.normalize();
    if (std::abs(n(2)) < thresh)
        return true;
    else
        return false;
}

float dist_to_plane(const Point &pt, const Plane &plane) {
    FT numerator = plane.a() * pt.x() + plane.b() * pt.y() +
                   plane.c() * pt.z() + plane.d();
    FT denominator_square =
        plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c();
    float dist = std::abs(CGAL::to_double(numerator)) /
                 std::sqrt(CGAL::to_double(denominator_square));
    return dist;
}



Eigen::Vector3d vector_of_two_point(
    const LCC_3::Point &p0,
    const LCC_3::Point &p1)
{
    return Eigen::Vector3d(CGAL::to_double(p1.x() - p0.x()), 
                           CGAL::to_double(p1.y() - p0.y()), 
                           CGAL::to_double(p1.z() - p0.z()));
}

bool on_same_line(
    const LCC_3 &lcc,
    Dart_const_descriptor &e1,
    Dart_const_descriptor &e2)
{
    LCC_3::Point p0 = lcc.point(e1), p1 = lcc.point(lcc.beta<1>(e1)),
                p2 = lcc.point(e2), p3 = lcc.point(lcc.beta<1>(e2));
    Eigen::Vector3d v1 = vector_of_two_point(p0, p1),
                    v2 = vector_of_two_point(p2, p3);
    v1.normalize();
    v2.normalize();

    if (std::abs(v1.dot(v2)) <= 1 - 1e-6)
        return false;
    else
    {
        if (lcc.info<0>(e1) != lcc.info<0>(e2))
        {
            Eigen::Vector3d v02 = vector_of_two_point(p0, p2);
            v02.normalize();
            if (std::abs(v1.dot(v02)) <= 1 - 1e-6)
                return false;
        }
        
        if (lcc.info<0>(e1) != lcc.info<0>(lcc.beta<2>(e2)))
        {
            Eigen::Vector3d v03 = vector_of_two_point(p0, p3);
            v03.normalize();
            if (std::abs(v1.dot(v03)) <= 1 - 1e-6)
                return false;
        }
    }
    return true;
}

std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

bool read_config(
    const char* path, 
    rj::Document& config_doc)
{
    FILE *config_fp = fopen(path, "rb");

    if (!config_fp)
    {
        std::cerr << "Error: unable to open config file\n";
        return false;
    }

    char config_readBuffer[65536];
    rj::FileReadStream config_is(config_fp, config_readBuffer,
                                 sizeof(config_readBuffer));
    config_doc.ParseStream(config_is);
    if (config_doc.HasParseError())
    {
        std::cerr << "Error: failed to parse JSON document\n";
        fclose(config_fp);
        return false;
    }
    fclose(config_fp);
    return true;
}

void make_combinations(
    const std::set<int> &numbers,
    std::back_insert_iterator<std::vector<std::pair<int, int>>> insertIt)
{
    for (auto it1 = numbers.begin(); it1 != numbers.end(); ++it1)
        for (auto it2 = std::next(it1); it2 != numbers.end(); ++it2)
            *(insertIt++) = *it1 < *it2 ? 
                            std::make_pair(*it1, *it2) 
                          : std::make_pair(*it2, *it1);
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> two_ends_of_edge(
    const LCC_3 &lcc,
    Dart_const_descriptor &edge)
{
    Eigen::Vector3d end1 (CGAL::to_double(lcc.point(edge).x()), 
                                    CGAL::to_double(lcc.point(edge).y()), 
                                    CGAL::to_double(lcc.point(edge).z())),
                    end2 (CGAL::to_double(lcc.point(lcc.beta<1>(edge)).x()),
                                    CGAL::to_double(lcc.point(lcc.beta<1>(edge)).y()),
                                    CGAL::to_double(lcc.point(lcc.beta<1>(edge)).z()));
    return std::make_pair(end1, end2);
}

Eigen::Vector3d inward_normal(
    const LCC_3 &lcc,
    const Dart_const_descriptor &face)
{
    Dart_const_descriptor this_edge = face;
    const Plane *plane1 = lcc.info<2>(this_edge).second,
                *plane2 = lcc.info<2>(lcc.beta<2>(this_edge)).second;
    Dart_const_descriptor next_edge = lcc.beta<1>(this_edge);
    while (on_same_line(lcc, next_edge, this_edge) && next_edge != face)
        next_edge = lcc.beta<1>(next_edge);
    if (next_edge == this_edge)
    {
        // draw
        std::cout << "error: no next edge\n";
        // #if VIS==1
        // LCC_3 draw_lcc = lcc;
        // size_type amark = draw_lcc.get_new_mark();
        // for (LCC_3::One_dart_per_cell_range<2>::iterator
        //     it=draw_lcc.one_dart_per_cell<2>().begin(),
        //     itend=draw_lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
        // {
        //     if (draw_lcc.info<2>(it).first == lcc.info<2>(face).first)
        //     {
        //         for (LCC_3::Dart_of_cell_range<2>::iterator
        //             itt=draw_lcc.darts_of_cell<2>(it).begin(),
        //             ittend=draw_lcc.darts_of_cell<2>(it).end(); itt!=ittend; ++itt)
        //             draw_lcc.mark(itt, amark);
        //     }
        // }
        // LCC_gs_options gso(amark);
        // CGAL::draw(draw_lcc, gso, "failed to calculate the inward normal!");
        // draw_lcc.free_mark(amark);
        // #endif
    }
    if (next_edge == this_edge) return Eigen::Vector3d(0, 0, 0);
    this_edge = lcc.beta<0>(next_edge);

    auto [end1_this_edge, end2_this_edge] = two_ends_of_edge(lcc, this_edge);
    auto [end1_next_edge, end2_next_edge] = two_ends_of_edge(lcc, next_edge);

    Eigen::Vector3d n_this_face = (end2_next_edge - end1_next_edge).cross(
                                    end1_this_edge - end2_this_edge);
    Eigen::Vector3d n(
        CGAL::to_double(plane1->a()), CGAL::to_double(plane1->b()), CGAL::to_double(plane1->c()));
    n.normalize();
    n = n.dot(n_this_face) > 0 ? n : -n;

    return n;
}

bool in_cell_set(
    const LCC_3 &lcc,
    const Dart_descriptor &cell,
    const std::vector<Dart_descriptor> &cell_set)
{
    for (const auto &c : cell_set)
    {
        if (lcc.info<3>(c).first == lcc.info<3>(cell).first)
            return true;
    }
    return false;
}

bool in_face_set(
    const LCC_3 &lcc,
    const Dart_descriptor &face,
    const std::set<Dart_descriptor> &face_set)
{
    for (const auto &f : face_set)
    {
        if (lcc.info<2>(f).first == lcc.info<2>(face).first)
            return true;
    }
    return false;
}

std::set<Volume_idx> nbh_cell3(
    const LCC_3 &lcc,
    const Volume_idx &vid)
{
    assert (vid >= 0 && vid < lcc.number_of_attributes<3>());
    std::set<Volume_idx> nbhs;
    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (lcc.info<3>(it).first == vid)            
        {
            #if DEBUG_CENTRAL==1
            std::cout << "finding nbh... \n";
            #endif

            int num_nbh = 0;
            for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator
                f=lcc.one_dart_per_incident_cell<2, 3>(it).begin();
                f!=lcc.one_dart_per_incident_cell<2, 3>(it).end(); ++f)
            {
                if (lcc.is_free<3>(f)) continue;
                Dart_const_descriptor oppo_f = lcc.beta<3>(f);
                Volume_idx this_v = lcc.info<3>(f).first;
                Volume_idx oppo_v = lcc.info<3>(oppo_f).first;
                bool oppo_in = lcc.info<3>(oppo_f).second;
                if (this_v != oppo_v && oppo_in)
                    nbhs.insert(oppo_v);
            }
            break;
        }
    }
    return nbhs;
}

int num_nbh_cell3(
    const LCC_3 &lcc,
    const Volume_idx &vid)
{
    std::set<Volume_idx> nbhs = nbh_cell3(lcc, vid);
    return nbhs.size();
}

float radian_from_u_to_v(
    const Eigen::Vector3d p0,
    const Eigen::Vector3d p1,
    const Eigen::Vector3d p2)
{
    Eigen::Vector3d u = p0 - p1, v = p2 - p1;
    Eigen::Vector3d cross = u.cross(v);
    float theta = std::atan2(cross.norm(), u.dot(v));

    Eigen::Vector3d n_up (0, 0, 1);
    if (cross.dot(n_up) < 0)
        theta = 2 * M_PI - theta;

    return theta;
}

float radian_from_u_to_v(
    const LCC_3::Point &p0,
    const LCC_3::Point &p1,
    const LCC_3::Point &p2,
    const Plane *plane)
{
    Point2 p0_2d = plane->to_2d(p0),
           p1_2d = plane->to_2d(p1),
           p2_2d = plane->to_2d(p2);
    Eigen::Vector3d p0_3d (CGAL::to_double(p0_2d.x()), CGAL::to_double(p0_2d.y()), 0),
                    p1_3d (CGAL::to_double(p1_2d.x()), CGAL::to_double(p1_2d.y()), 0),
                    p2_3d (CGAL::to_double(p2_2d.x()), CGAL::to_double(p2_2d.y()), 0);
    return radian_from_u_to_v(p0_3d, p1_3d, p2_3d);
}


bool counterclockwise(
    const Eigen::Vector3d p0,
    const Eigen::Vector3d p1,
    const Eigen::Vector3d p2)
{
    Eigen::Vector3d u = p0 - p1, v = p2 - p1;
    u.normalize();
    v.normalize();
    Eigen::Vector3d cross = u.cross(v);
    if (cross.norm() < 1e-3)
        return true;
    return cross.z() < 0; 
}

bool is_concave(
    const LCC_3 &lcc,
    const Dart_const_descriptor &p0,
    const Dart_const_descriptor &p1,
    const Dart_const_descriptor &p2,
    const Plane *plane)
{
    LCC_3::Point p0_3d = lcc.point(p0),
                 p1_3d = lcc.point(p1),
                 p2_3d = lcc.point(p2);
    Point2 p0_2d = plane->to_2d(p0_3d),
            p1_2d = plane->to_2d(p1_3d),
            p2_2d = plane->to_2d(p2_3d);  
    Eigen::Vector3d p0_3d_v (CGAL::to_double(p0_2d.x()), CGAL::to_double(p0_2d.y()), 0),
                    p1_3d_v (CGAL::to_double(p1_2d.x()), CGAL::to_double(p1_2d.y()), 0),
                    p2_3d_v (CGAL::to_double(p2_2d.x()), CGAL::to_double(p2_2d.y()), 0);
    Eigen::Vector3d v10 = p0_3d_v - p1_3d_v,
                    v12 = p2_3d_v - p1_3d_v;
    v10.normalize();
    v12.normalize();
    float norm = v10.cross(v12).norm();
    if (norm < 1e-6)
        return false;
    return v10.cross(v12).z() > 0;
}

bool counterclockwise(
    const std::vector<Point2> &pts)
{
    FT sum = 0;
    for (int i = 0; i < pts.size(); i++)
        sum += (pts[(i+1)%pts.size()].x() - pts[i].x()) * (pts[(i+1)%pts.size()].y() + pts[i].y());
    return sum < 0;
}

bool counterclockwise(
    const LCC_3 &lcc,
    const std::vector<Dart_const_descriptor> &cycle,
    const Plane *plane)
{
    std::vector<Point2> pts;
    for (const auto &v : cycle)
        pts.push_back(plane->to_2d(lcc.point(v)));
    return counterclockwise(pts);
}

bool counterclockwise(
    const LCC_3 &lcc, 
    const std::vector<int> &cycle,
    const std::map<std::pair<int, int>, Dart_descriptor> &edge_darts,
    const Plane *plane)
{
    std::vector<Dart_const_descriptor> cycle_darts;
    for (int i = 0; i < cycle.size(); i++)
    {
        int end1 = cycle[i], end2 = cycle[(i+1)%cycle.size()];
        std::pair<int, int> edge = end1 < end2 ? std::make_pair(end1, end2) : std::make_pair(end2, end1);
        Dart_const_descriptor edge_dart = end1 < end2 ? edge_darts.at(edge) : lcc.beta<2>(edge_darts.at(edge));
        cycle_darts.push_back(edge_dart);
    }

    return counterclockwise(lcc, cycle_darts, plane);
}

int break_concave_loop(
    LCC_3 &lcc,
    const std::vector<Dart_descriptor> &cycle,
    const Plane *plane,
    std::vector<std::pair<Dart_descriptor, Dart_descriptor>> &breakers)
{
    std::vector<Point2> pts;
    for (size_t i = 0; i < cycle.size(); i++)
    {
        Point2 p = plane->to_2d(lcc.point(cycle[i]));
        pts.push_back(p);
    }
    Partition_traits_2 traits(CGAL::make_property_map(pts));
    size_t num_pts = pts.size();

    Polygon2 plg;
    for (size_t i = 0; i < num_pts; i++)
        plg.push_back(i);
    if (!CGAL::is_simple_2(pts.begin(), pts.end())) 
    {
        #if DEBUG_CUT == 1
        std::cout << "not a simple polygon:";
        for (auto dart : cycle)
            std::cout << lcc.info<0>(dart) << "-";
        std::cout << "\n";
        #endif
        return 0;
    }
    
    std::list<Polygon2> partitions;
    CGAL::approx_convex_partition_2(plg.vertices_begin(),
                                plg.vertices_end(),
                                std::back_inserter(partitions), traits);

    std::set<std::pair<int, int>> break_edges;
    #if DEBUG_CUT == 1
    std::cout << "num of partitions: " << partitions.size() << "\n";
    #endif
    for(const auto par : partitions)
    {
        for (size_t i = 0; i < par.size(); i++)
        {
            int e1 = par[i], e2 = par[(i+1)%par.size()];
            if (std::abs(e2 - e1) != 1 && std::abs(e2 - e1) != num_pts - 1)
            {
                #if DEBUG_CUT == 1
                std::cout << "find a breaker!\n";
                #endif

                if (e1 < e2)
                    break_edges.insert(std::make_pair(e1, e2));
                else
                    break_edges.insert(std::make_pair(e2, e1));
                break;
            }
        }
    }
    
    for (const auto &[e1, e2] : break_edges)
    {
        #if DEBUG_CUT == 1
        std::cout << "breaker: " << lcc.info<0>(cycle[e1]) << "," << lcc.info<2>(cycle[e1]).first
                  << "-" 
                  << lcc.info<0>(cycle[e2]) << "," << lcc.info<2>(cycle[e2]).first << "\n";
        #endif
        breakers.push_back(std::make_pair(cycle[e1], cycle[e2]));
    }
    return breakers.size();
}

int get_cycles_from_arrangement(
    LCC_3 &lcc,
    const std::map<std::pair<int, int>, Dart_descriptor> &conn_inters,
    const Plane* plane,
    std::vector<std::vector<Dart_descriptor>> &cycles,
    Arrangement_2 &arr)
{
    // create arrangement
    std::vector<CM_Segment_2> segs;
    std::map<int, Point2> pts;
    for (auto const &[inter_pair, edge] : conn_inters)
    {
        int end1 = inter_pair.first, end2 = inter_pair.second;
        LCC_3::Point p0 = lcc.point(edge),
                     p1 = lcc.point(lcc.beta<1>(edge));
        Point2 p0_2d = plane->to_2d(p0),
               p1_2d = plane->to_2d(p1);
        pts[end1] = p0_2d;
        pts[end2] = p1_2d;
    }

    for (auto const &[inter_pair, edge] : conn_inters)
    {
        Dart_descriptor dart = edge, oppo_dart = lcc.beta<2>(edge);

        if (pts[inter_pair.first].x() < pts[inter_pair.second].x())
        {
            segs.push_back(CM_Segment_2(
                Segment_2(pts[inter_pair.first], pts[inter_pair.second]), dart));
        }
        else if (pts[inter_pair.first].x() == pts[inter_pair.second].x())
        {
            if (pts[inter_pair.first].y() < pts[inter_pair.second].y())
                segs.push_back(CM_Segment_2(
                    Segment_2(pts[inter_pair.first], pts[inter_pair.second]), dart));
            else if (pts[inter_pair.first].y() > pts[inter_pair.second].y())
                segs.push_back(CM_Segment_2(
                    Segment_2(pts[inter_pair.second], pts[inter_pair.first]), oppo_dart));
            else
                continue; // degenerated segment
        }
        else
            segs.push_back(CM_Segment_2(
                Segment_2(pts[inter_pair.second], pts[inter_pair.first]), oppo_dart));
    }
    CGAL::insert<Traits_2>(arr, segs.begin(), segs.end());
    // CGAL::insert(arr, segs.begin(), segs.end());

    //get bounded ccb faces
    Arrangement_2::Face_const_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {
        if (fit->is_unbounded()) continue;
        if (fit->number_of_holes() > 0) continue;
        std::vector<Dart_descriptor> cycle;
        Arrangement_2::Ccb_halfedge_const_circulator start, curr;
        start = curr = fit->outer_ccb();

        do{ 
            Arrangement_2::Halfedge_const_handle he = curr;
            if (curr->curve().data().size() != 1)
                std::cout << "curr->curve().data().size() ==" << curr->curve().data().size() << "\n";
            Dart_descriptor dart = curr->curve().data().front(); 
            // dart is for (source.x < target.x) || (source.x == target.x && source.y < target.y)
            // if source and target don't satisfy this, should inverse the dart
            if (he->source()->point().x() < he->target()->point().x())
                cycle.push_back(dart);
            else if (he->source()->point().x() == he->target()->point().x())
            {
                if (he->source()->point().y() < he->target()->point().y())
                    cycle.push_back(dart);
                else
                    cycle.push_back(lcc.beta<2>(dart));
            }
            else
                cycle.push_back(lcc.beta<2>(dart));
        }while(++curr != start);

        if (cycle.size() > 0)
            cycles.push_back(cycle);
    }
    return cycles.size();
}

 
int trace_loops(
    LCC_3 &lcc, 
    const std::map<std::pair<int, int>, Dart_descriptor> &conn_inters,
    const Plane* cut_plane,
    std::vector<std::pair<std::vector<Dart_descriptor>, std::vector<std::pair<Dart_descriptor, Dart_descriptor>>>> &filtered_cycles,
    Arrangement_2 &arr)
{   
    std::vector<std::vector<Dart_descriptor>> cycles;
    get_cycles_from_arrangement(lcc, conn_inters, cut_plane, cycles, arr);
    #if DEBUG_CUT==1
    std::cout << "arrangement #cycles=" << cycles.size() << "\n";
    #endif

    for (auto cycle : cycles)
    {
        bool on_plane_face = true, concave = false;
        std::set<int> start_adj;
        for (int i = 0; i < cycle.size(); i++)
        {
            Dart_descriptor this_dart = cycle[i],
                            oppo_dart = lcc.beta<2>(cycle[i]),
                            next_part = cycle[(i+1)%cycle.size()],
                            prev_part = cycle[(i-1+cycle.size())%cycle.size()];
            
            // if (lcc.info<2>(this_dart).second != cut_plane && lcc.info<2>(oppo_dart).second != cut_plane)
            //     on_plane_face = false;
            if (i == 0)
                start_adj = {lcc.info<2>(this_dart).first, lcc.info<2>(oppo_dart).first};
            if (i > 0)
                if (start_adj.find(lcc.info<2>(this_dart).first) == start_adj.end()
                    && start_adj.find(lcc.info<2>(oppo_dart).first) == start_adj.end())
                    on_plane_face = false;

            if (is_concave(lcc, prev_part, this_dart, next_part, cut_plane))
                concave = true;
        }

        if (on_plane_face) 
        {
            #if DEBUG_CUT == 1
            std::cout << "[CUT] filtered out: on plane=" << on_plane_face << ", concave=" << concave << "\n";
            for (auto d : cycle)
                std::cout << lcc.info<0>(d) << "(" << lcc.info<3>(d).first << ") ";
            std::cout << "\n";
            #endif
    
            continue;
        }
        std::vector<std::pair<Dart_descriptor, Dart_descriptor>> concave_breakers;
        if (concave)
        {
            int num_breaker = break_concave_loop(lcc, cycle, cut_plane, concave_breakers);
            #if DEBUG_CUT == 1
            std::cout << "concave, found #breaker=" << num_breaker << "\n";
            #endif
            if (num_breaker == 0) continue;
        }

        filtered_cycles.push_back(std::make_pair(cycle, concave_breakers));
    }

    #if DEBUG_CUT == 1
    std::cout << "num of filtered cycles: " << filtered_cycles.size() << "\n";
    for (int i = 0; i < filtered_cycles.size(); i++)
    {
        std::cout << "[CUT] cycle " << i << ": ";
        for (auto d : filtered_cycles[i].first)
            std::cout << lcc.info<0>(d) << " ";
        std::cout << "\n";
    }
    #endif
    
    return filtered_cycles.size();
}

int cut_cell3_by_plane(
    LCC_3 &lcc,
    Dart_descriptor &cell_dart,
    const Plane *plane,
    std::vector<Dart_descriptor> &pos_cells,
    std::vector<Dart_descriptor> &neg_cells,
    bool convex=false)
{
    Eigen::Vector3d plane_normal(
        CGAL::to_double(plane->a()), CGAL::to_double(plane->b()), CGAL::to_double(plane->c()));
    // 1. insert intersection points (0-cells)

    // 1.1 get all edges of the cell
    std::unordered_map<int, Dart_descriptor> darts_of_cell_edges;
    for (LCC_3::One_dart_per_incident_cell_range<1, 3>::iterator
        it=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).begin(),
        itend=lcc.one_dart_per_incident_cell<1, 3>(cell_dart).end(); it!=itend; ++it)
        darts_of_cell_edges[lcc.info<1>(it)] = it;
    
    #if DEBUG_CUT == 1
    std::cout << "#edges=" << darts_of_cell_edges.size() << " in cell3-" << lcc.info<3>(cell_dart).first << "\n";
    #endif

    std::map<int, Dart_descriptor> inters;
    std::map<int, const Plane*> adj_face_darts;
    std::map<int, std::set<int>> vertices_on_faces; //loop the 2-cells of each vertices
    std::map<std::pair<int, int>, Dart_descriptor> conn_inters;
    for (const auto [eid, e] : darts_of_cell_edges)
    {
        // 1.2 Get the segment of the edge
        Dart_descriptor vtx0 = e;
        Dart_descriptor vtx1 = lcc.beta<2>(e);
        LCC_3::Point p0 = lcc.point(vtx0);
        LCC_3::Point p1 = lcc.point(vtx1);
        Segment seg(Point(p0.x(), p0.y(), p0.z()), Point(p1.x(), p1.y(), p1.z()));
        
        // 1.3. Get the intersection point
        if (CGAL::do_intersect(seg, *plane))
        {
            std::pair<int, int> adj_faces = 
                    std::make_pair(lcc.info<2>(vtx0).first, lcc.info<2>(vtx1).first);
            adj_face_darts[adj_faces.first] = lcc.info<2>(e).second;
            adj_face_darts[adj_faces.second] = lcc.info<2>(lcc.beta<2>(e)).second;

            CGAL::Object obj = CGAL::intersection(seg, *plane);
            const Point* p = CGAL::object_cast<Point>(&obj);
            const Segment* s = CGAL::object_cast<Segment>(&obj);
            // std::vector<std::pair<int, Dart_descriptor>> edge_inters;
            std::unordered_map<int, Dart_descriptor> edge_inters;

            if (p)
            {
                if (*p == seg.source())
                    inters[lcc.info<0>(vtx0)] = vtx0;
                else if (*p == seg.target())
                    inters[lcc.info<0>(vtx1)] = vtx1;
                else
                {
                    // creat a new vertex on the edge
                    Dart_descriptor d = lcc.insert_point_in_cell<1>(e, LCC_3::Point(p->x(), p->y(), p->z())); // a dart originated from the new vertex                    
                    lcc.info<0>(d) = lcc.number_of_attributes<0>() - 1;
                    assert(lcc.info<3>(d) == lcc.info<3>(cell_dart));
                    inters[lcc.info<0>(d)] = d;
                    #if DEBUG_CUT == 1
                    std::cout << "[CUT] inserted a new vertex\n";
                    #endif
                }
            }
            else if (s)
            {
                std::pair<int, int> inter_pair = lcc.info<0>(vtx0) < lcc.info<0>(vtx1) ? 
                                                std::make_pair(lcc.info<0>(vtx0), lcc.info<0>(vtx1)) 
                                              : std::make_pair(lcc.info<0>(vtx1), lcc.info<0>(vtx0));
                conn_inters[inter_pair] = lcc.info<0>(vtx0) < lcc.info<0>(vtx1) ? vtx0: vtx1;
                inters[lcc.info<0>(vtx0)] = vtx0;
                inters[lcc.info<0>(vtx1)] = vtx1;
            }
            else
            {
                std::cerr << "\t\tERROR: intersection is not a point or a segment\n";
                return EXIT_FAILURE;
            }
            // add the intersection of this edge to inters and vertices_on_faces
        }
    }
    std::map<int, const Plane*> face_plane;
    for (const auto [vtx_id, vtx_dart] : inters)
    {
        #if DEBUG_CUT == 1
        std::cout << "[CUT] vtx_id=" << vtx_id << " incident cell2: ";
        #endif

        for (LCC_3::One_dart_per_incident_cell_range<2, 0>::iterator
            it=lcc.one_dart_per_incident_cell<2, 0>(vtx_dart).begin(),
            itend=lcc.one_dart_per_incident_cell<2, 0>(vtx_dart).end(); it!=itend; ++it)
        {
            if (lcc.info<3>(it) != lcc.info<3>(cell_dart)) continue;
            #if DEBUG_CUT == 1
            std::cout << lcc.info<2>(it).first << " ";
            #endif
            vertices_on_faces[lcc.info<2>(it).first].insert(vtx_id);
            face_plane[lcc.info<2>(it).first] = lcc.info<2>(it).second;
        }
        #if DEBUG_CUT == 1
        std::cout << "\n";
        #endif
    }

    #if DEBUG_CUT == 1
    std::cout << "#inters=" << inters.size() 
              << ", #intersected faces=" << vertices_on_faces.size() << "\n";
    if (!convex)
    {
        // #if VIS_CONCAVE_CUT == 1 && VIS == 1
        // CGAL::draw(lcc, "inserted vertices");
        // #endif
    }
    #endif

    // 2. connect edges (1-cells) between intersection points on the same face
    for (const auto [fid, f_inters] : vertices_on_faces)
    {
        if (face_plane[fid] == plane) continue;
        #if DEBUG_CUT == 1
        std::cout << "face " << fid << " inters: ";
        for (auto inter : f_inters)
            std::cout << inter << "-";
        std::cout << "\n";
        #endif
        // faces created in the reconstruction are all convex, simplifying the connecting between intersections        
        std::vector<std::pair<int, int>> inter_combo;
        make_combinations(f_inters, std::back_inserter(inter_combo));
        bool exist_edges = false;
        for (const auto inter: inter_combo)
        {
            if (conn_inters.find(inter) != conn_inters.end())
            {
                exist_edges = true;
                #if DEBUG_CUT == 1
                std::cout << "existing an connected edge on the face, skip it \n";
                #endif
                break;
            }
        }
        if (exist_edges) continue;

        for (const auto [inter1, inter2] : inter_combo)
        {
            #if DEBUG_CUT == 1
            std::cout << "[CUT] edge: " << inter1 << "(" <<  lcc.info<2>(inters[inter1]).first << ")-"
                                        << inter2 << "(" <<  lcc.info<2>(inters[inter2]).first << ")\n";
            #endif

            Dart_descriptor dart1, dart2;
            for (LCC_3::One_dart_per_incident_cell_range<2, 0>::iterator
                it=lcc.one_dart_per_incident_cell<2, 0>(inters[inter1]).begin(),
                itend=lcc.one_dart_per_incident_cell<2, 0>(inters[inter1]).end(); it!=itend; ++it)
            {
                if (lcc.info<2>(it).first == fid)
                {
                    dart1 = it;
                    break;
                }
            }

            for (LCC_3::One_dart_per_incident_cell_range<2, 0>::iterator
                it=lcc.one_dart_per_incident_cell<2, 0>(inters[inter2]).begin(),
                itend=lcc.one_dart_per_incident_cell<2, 0>(inters[inter2]).end(); it!=itend; ++it)
            {
                if (lcc.info<2>(it).first == fid)
                {
                    dart2 = it;
                    break;
                }
            }
            if (dart1 == nullptr || dart2 == nullptr)
            {
                std::cerr << "\t\tERROR: edge not found\n"; // why??????
                // size_type amark = lcc.get_new_mark();
                // for (LCC_3::Dart_of_cell_range<2>::const_iterator
                //     it=lcc.darts_of_cell<2>(e_dart).begin(),
                //     itend=lcc.darts_of_cell<2>(e_dart).end(); it!=itend; ++it)
                //     lcc.mark(it, amark);

                continue;
            }
            #if DEBUG_CUT == 1
            std::cout << "[CUT] dart1=" << lcc.info<2>(dart1).first << "-dart2=" << lcc.info<2>(dart2).first << "\n";
            #endif
            
            // skip connections between points on a same line (intersection between a same pair of planes)
            // if (lcc.info<2>(inters[inter1]).first != lcc.info<2>(inters[inter2]).first)
            // {
            //     if (lcc.info<1>(inters[inter1]) == lcc.info<1>(inters[inter2]))
            //     {
            //         std::cout << "on the same edge\n";
            //         conn_inters[std::make_pair(inter1, inter2)] = inters[inter1];
            //         continue;
            //     }
            // }
            
            if (lcc.beta<1>(dart1) == dart2)
            {
                #if DEBUG_CUT == 1
                std::cout << "[CUT] on the same edge\n";
                #endif 

                conn_inters[std::make_pair(inter1, inter2)] = dart1;
                continue;
            }
            if (lcc.beta<0>(dart1) == dart2)
            {
                #if DEBUG_CUT == 1
                std::cout << "[CUT] on the same edge\n";
                #endif

                conn_inters[std::make_pair(inter1, inter2)] = lcc.beta<0, 2>(dart1);
                continue;
            }
            if (lcc.info<2>(lcc.beta<2>(dart1)).second == lcc.info<2>(lcc.beta<2>(dart2)).second)
            {
                // 在一个平面上共beta2面也不一定在same line上，有种可能是这三个面都共面
                Dart_const_descriptor e1 = dart1, e2 = dart2;
                if (on_same_line(lcc, e1, e2))
                {
                    #if DEBUG_CUT == 1
                    std::cout << "[CUT] intersecion points on the same line but not on the same edge\n";
                    #endif
                    continue;
                }
            }
            
            if (lcc.is_insertable_cell_1_in_cell_2(dart1, dart2)) // check if dart1 belong to beta<1> orbit of dart2
            {
                Dart_descriptor new_edge = lcc.insert_cell_1_in_cell_2(dart1, dart2);
                lcc.set_attribute<1>(new_edge, lcc.create_attribute<1>());
                lcc.info<1>(new_edge) = lcc.number_of_attributes<1>() - 1;
                conn_inters[std::make_pair(inter1, inter2)] = lcc.beta<2>(new_edge);
                #if DEBUG_CUT == 1
                std::cout << "[CUT] inserted an edge\n";
                #endif
            }
            else
            {
                #if DEBUG_CUT == 1
                std::cerr << "[CUT] ERROR: edge insertion failed\n";
                std::cerr << "[CUT] reverse? " << lcc.is_insertable_cell_1_in_cell_2(dart2, dart1) << "\n";
                #endif
            }
                    
        }
    }
    #if DEBUG_CUT == 1
    std::cout << "[CUT] #conn_inters=" << conn_inters.size() << "\n";
    if (!convex)
    {
        // #if VIS_CONCAVE_CUT == 1 && VIS == 1
        // CGAL::draw(lcc, "inserted edges");
        // #endif
    }
    #endif
    // 3. detect the cut loops
    // std::set<Edge, EdgeCmp> edges;
    // for (const auto [inter_pair, edge] : conn_inters)
    //     edges.insert(
    //         {inter_pair.first, inter_pair.second, length_of_edge(lcc, edge)});
    
    // Graph g = Graph(edges);
    // std::vector<std::vector<int>> cycles;
    // g.detect_all_cycle_bases(cycles);
    // #if DEBUG_CUT == 1
    // std::cout << "#cycles=" << cycles.size() << "\n";
    // #endif
    // if (cycles.size() == 0)
    // {
    //     #if DEBUG_CUT == 1
    //     std::cerr << "\t\tERROR: no cut loops detected\n";
    //     #endif
    //     return EXIT_FAILURE;
    // }
    // if (convex)
    // {
    //     if (cycles.size() > 1)
    //     {
    //         #if DEBUG_CUT == 1
    //         std::cerr << "\t\tERROR: more than one cut loops detected in a convex cell, select the largest one only\n";
    //         #endif
    //         for (int i = 1; i < cycles.size(); i++)
    //             cycles.pop_back(); 
    //     }
    // }
    std::vector<std::pair<std::vector<Dart_descriptor>, std::vector<std::pair<Dart_descriptor, Dart_descriptor>>>> cycles;
    Arrangement_2 arr;
    int num_cycle = trace_loops(lcc, conn_inters, plane, cycles, arr);

    #if DEBUG_CUT == 1
    std::cout << "num_cycle: " << num_cycle << "\n";
    #endif
    
    if (cycles.size() == 0 && inters.size() > 0 && conn_inters.size() > 0)
    {
        #if VIS_CONCAVE_CUT == 1 && VIS == 1
        CGAL::draw(arr, "Arrangement of no loops detected");
        size_type amark = lcc.get_new_mark();
        for (LCC_3::One_dart_per_cell_range<2>::iterator
            it=lcc.one_dart_per_cell<2>().begin(),
            itend=lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
        {
            if (lcc.info<2>(it).second == plane)
            {
                for (LCC_3::Dart_of_cell_range<2>::iterator
                    itt=lcc.darts_of_cell<2>(it).begin(),
                    ittend=lcc.darts_of_cell<2>(it).end(); itt!=ittend; ++itt)
                    lcc.mark(itt, amark);
            }
        }
        LCC_gs_options gso(amark);
        CGAL::draw(lcc, gso, "no cut loops detected");
        lcc.free_mark(amark);
        #endif
    }

    // 4. insert faces (2-cells)
    // reminder: attach the plane with the newly-inserted faces
    for (int i = 0; i < cycles.size(); i++)
    {
        // converting the cycle into ordered darts
        std::vector<Dart_descriptor> cycle_darts = cycles[i].first;
        std::vector<std::pair<Dart_descriptor, Dart_descriptor>> ccv_breakers = cycles[i].second;

        if (lcc.is_insertable_cell_2_in_cell_3(cycle_darts.begin(), cycle_darts.end()))
        {
            Dart_descriptor new_face = 
                    lcc.insert_cell_2_in_cell_3(cycle_darts.begin(), cycle_darts.end());
            lcc.set_attribute<2>(new_face, lcc.create_attribute<2>());
            lcc.info<2>(new_face) = std::make_pair(lcc.number_of_attributes<2>() - 1, plane);
            Dart_descriptor oppo_new_face = lcc.beta<3>(new_face);
            assert(lcc.info<2>(new_face) == lcc.info<2>(oppo_new_face));

            if (lcc.info<3>(new_face) == lcc.info<3>(oppo_new_face))
            {
                #if DEBUG_CUT == 1
                std::cout << "[CUT] lcc.info<3>(new_face) == lcc.info<3>(oppo_new_face)\n";
                #endif
                if (convex)
                {
                    lcc.remove_cell<2>(new_face);
                    return EXIT_FAILURE;
                }
                else
                    continue; // when cutting a ring, two cycles are needed for spearate the two parts
            }
            #if DEBUG_CUT == 1
            std::cout << "[CUT] inserted a face\n";
            #endif
            Eigen::Vector3d n1 = inward_normal(lcc, new_face);
            if (n1.norm() < 1e-3) continue;
             Eigen::Vector3d n2 = -n1;
            if (!in_cell_set(lcc, new_face, pos_cells) && !in_cell_set(lcc, new_face, neg_cells))
            {
                if (n1.dot(plane_normal) > 0)
                    pos_cells.push_back(new_face);
                else
                    neg_cells.push_back(new_face);
            }

            if (!in_cell_set(lcc, oppo_new_face, pos_cells) && !in_cell_set(lcc, oppo_new_face, neg_cells))
            {
                if (n2.dot(plane_normal) > 0)
                    pos_cells.push_back(oppo_new_face);
                else
                    neg_cells.push_back(oppo_new_face);
            }

            if (ccv_breakers.size() > 0)
            {
                #if DEBUG_CUT == 1
                std::cout << "inserted a concave face and let's break it into convex ones!\n";
                #endif
                for (auto [e1, e2] : ccv_breakers)
                {
                    Dart_descriptor insert_d1, insert_d2;
                    for (LCC_3::One_dart_per_incident_cell_range<2, 0>::iterator
                        it1=lcc.one_dart_per_incident_cell<2, 0>(e1).begin(),
                        it1end=lcc.one_dart_per_incident_cell<2, 0>(e1).end(); it1!=it1end; ++it1)
                    {
                        for (LCC_3::One_dart_per_incident_cell_range<2, 0>::iterator
                            it2=lcc.one_dart_per_incident_cell<2, 0>(e2).begin(),
                            it2end=lcc.one_dart_per_incident_cell<2, 0>(e2).end(); it2!=it2end; ++it2)
                        {
                            if (lcc.info<2>(it1).first == lcc.info<2>(it2).first)
                            {
                                insert_d1 = it1;
                                insert_d2 = it2;
                                break;
                            }
                        }
                    }

                    // Dart_descriptor insert_d1 = lcc.beta<2, 1>(e1), insert_d2 = lcc.beta<2, 1>(e2);
                    
                    // if (lcc.info<2>(insert_d2).first != lcc.info<2>(insert_d1).first)
                    // {
                    //     while (lcc.info<2>(insert_d2).first != lcc.info<2>(insert_d1).first && insert_d2 != e2)
                    //         insert_d2 = lcc.beta<2, 1>(insert_d2);
                    // }
                    
                    #if DEBUG_CUT == 1
                    std::cout << "e1: " << lcc.info<0>(insert_d1) << "(" << lcc.info<2>(insert_d1).first << "," << lcc.info<3>(insert_d1).first << ")"
                             << ", e2: " << lcc.info<0>(insert_d2) << "(" << lcc.info<2>(insert_d2).first << "," << lcc.info<3>(insert_d2).first << ")\n";
                    #endif

                   if (lcc.is_insertable_cell_1_in_cell_2(insert_d1, insert_d2))
                    {
                        Dart_descriptor new_edge = lcc.insert_cell_1_in_cell_2(insert_d1, insert_d2);
                        lcc.set_attribute<1>(new_edge, lcc.create_attribute<1>());
                        lcc.info<1>(new_edge) = lcc.number_of_attributes<1>() - 1;
                        #if DEBUG_CUT == 1
                        std::cout << "[CUT] inserted an edge for breaking concave face\n";
                        #endif
                    }
                    else
                    {
                        #if DEBUG_CUT == 1
                        std::cerr << "[CUT] ERROR: edge insertion for concave breaking failed\n";
                        #endif
                    }
                }
            }
        }
        else
        {   
            #if DEBUG_CUT == 1
            std::cerr << "[CUT] ERROR: face insertion failed\n";
            #endif
        }
            
    }
    #if DEBUG_CUT == 1
    std::cout << "[CUT] #pos_cells=" << pos_cells.size() << ", #neg_cells=" << neg_cells.size() << "\n";
    if (!convex)
    {
        if (pos_cells.size() == 0 || neg_cells.size() == 0)
        {
            #if VIS_CONCAVE_CUT == 1 && VIS == 1
            CGAL::draw(lcc, "inserted faces failed");
            #endif
        }
    }
    #endif

    if (pos_cells.size() > 0 && neg_cells.size() > 0)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

float inner_dihedral_angle(
    const LCC_3 &lcc,
    Dart_const_descriptor &edge_of_the_dart, 
    Dart_const_descriptor &cell_dart,
    bool &flip1, bool &flip2)
{

    Volume_idx vid = lcc.info<3>(cell_dart).first;
    assert (lcc.info<3>(edge_of_the_dart).first == vid);
    Dart_const_descriptor oppo_edge = lcc.beta<2>(edge_of_the_dart);

    const Plane *plane1 = lcc.info<2>(edge_of_the_dart).second,
                *plane2 = lcc.info<2>(oppo_edge).second;
    Eigen::Vector3d n1 (CGAL::to_double(plane1->a()), CGAL::to_double(plane1->b()), CGAL::to_double(plane1->c())),
                    n2 (CGAL::to_double(plane2->a()), CGAL::to_double(plane2->b()), CGAL::to_double(plane2->c()));
    n1.normalize();
    n2.normalize();

    if (plane1 == plane2) return M_PI;
    Eigen::Vector3d n1_inward = inward_normal(lcc, edge_of_the_dart),
                    n2_inward = inward_normal(lcc, oppo_edge);
    if (n1_inward.norm() < 1e-3 || n2_inward.norm() < 1e-3)
    {
        std::cerr << "\t\tERROR: inward normal is zero\n";
        return -1.0;
    }
    n1 = n1.dot(n1_inward) > 0 ? n1 : -n1;
    n2 = n2.dot(n2_inward) > 0 ? n2 : -n2;
    flip1 = n1.dot(n1_inward) < 0;
    flip2 = n2.dot(n2_inward) < 0;
    #if DEBUG_CONCAVE_EDGE == 1
    std::cout << "n1: " << n1.transpose() << ", n2: " << n2.transpose() << "\n";
    #endif

    float theta = std::acos(n1.dot(n2));

    Eigen::Vector3d end1_edge_face1 (CGAL::to_double(lcc.point(edge_of_the_dart).x()),
                                     CGAL::to_double(lcc.point(edge_of_the_dart).y()),
                                     CGAL::to_double(lcc.point(edge_of_the_dart).z())),
                    end2_edge_face1 (CGAL::to_double(lcc.point(lcc.beta<1>(edge_of_the_dart)).x()),
                                     CGAL::to_double(lcc.point(lcc.beta<2>(edge_of_the_dart)).y()),
                                     CGAL::to_double(lcc.point(lcc.beta<2>(edge_of_the_dart)).z())),
                    end1_edge_face2 (CGAL::to_double(lcc.point(oppo_edge).x()),
                                     CGAL::to_double(lcc.point(oppo_edge).y()),
                                     CGAL::to_double(lcc.point(oppo_edge).z())),
                    end2_edge_face2 (CGAL::to_double(lcc.point(lcc.beta<1>(oppo_edge)).x()),
                                     CGAL::to_double(lcc.point(lcc.beta<1>(oppo_edge)).y()),
                                     CGAL::to_double(lcc.point(lcc.beta<1>(oppo_edge)).z()));

    Eigen::Vector3d f1 = n1.cross(end2_edge_face1 - end1_edge_face1),
                    f2 = n2.cross(end2_edge_face2 - end1_edge_face2);
    f1.normalize();
    f2.normalize();
    
    #if DEBUG_CONCAVE_EDGE == 1
    std::cout << "f1: " << f1.transpose() << ", f2: " << f2.transpose() << "\n";
    #endif

    bool concave = (f1 + f2).dot(n1 + n2) < 0;
    float alpha = concave ? M_PI + theta : M_PI - theta;

    return alpha;
}


void extract_triangle_mesh(
    const LCC_3 &lcc, 
    const Volume_idx &vid,
    std::vector<Point> &pts, 
    std::vector<std::vector<size_t>> &triangles,
    std::function<Point(Point)> transform = nullptr) 
{
    assert (vid >= 0 && vid < lcc.number_of_attributes<3>());
    Dart_const_descriptor cell_dart;
    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (lcc.info<3>(it).first == vid)            
        {
            cell_dart = it;
            break;
        }
    } 

    pts.clear();
    triangles.clear();
    std::vector<int> all_lcc_pts;
    std::vector<int>::iterator lcc_pt_it;
    std::vector<int> all_vids;

    for (LCC_3::One_dart_per_incident_cell_range<2, 3>::const_iterator
            it=lcc.one_dart_per_incident_cell<2, 3>(cell_dart).begin(),
            itend=lcc.one_dart_per_incident_cell<2, 3>(cell_dart).end(); it!=itend; ++it)
    {
        // surface face
        std::vector<size_t> face;
        Dart_const_descriptor d = it;
        do
        {
            LCC_3::Point lcc_pt = lcc.point(d);
            int vtx_idx = lcc.info<0>(d);
            lcc_pt_it = std::find(all_lcc_pts.begin(), all_lcc_pts.end(), vtx_idx);
            if (lcc_pt_it == all_lcc_pts.end())
            {
                Point pt(lcc_pt.x(), lcc_pt.y(), lcc_pt.z());
                if (transform != nullptr)
                    pt = transform(pt);
                size_t vid = pts.size();
                all_lcc_pts.push_back(vtx_idx);
                all_vids.push_back(vid);
                pts.push_back(pt);
                face.push_back(vid);
            }
            else
            {
                int idx = std::distance(all_lcc_pts.begin(), lcc_pt_it);
                face.push_back(all_vids[idx]);
            }
            d = lcc.beta<1>(d);
        } while (d != it);
        // faces.push_back(face);
        
        for (int i = 1; i < face.size() - 1; i++)
            triangles.push_back({face[0], face[i], face[i + 1]});
    }
    // created through mesh
    // std::cout << "creating surface mesh ... \n";
    // // std::cout << "pts size: " << pts.size() << " faces size: " << faces.size() << std::endl;
    // PMP::orient_polygon_soup(pts, faces);
    // std::cout << "orient done\n";
    // PMP::repair_polygon_soup(pts, faces); // required for robust ray intersection
    // std::cout << "repaired done\n";
    // PMP::polygon_soup_to_polygon_mesh(pts, faces, mesh);
    // std::cout << "polygon_soup_to_polygon_mesh done\n";
    // PMP::stitch_borders(mesh);
    // std::cout << "stitch_borders done\n";
    // PMP::triangulate_faces(mesh);
    // std::cout << "triangulate_faces done\n";
    // assert(CGAL::is_closed(mesh));
    // std::cout << "closed mesh!\n";
    // PMP::orient_to_bound_a_volume(mesh);
    // std::cout << "orient_to_bound_a_volume done\n";
    // CGAL::draw(mesh);
    // // ----- created through surface mesh
}
 

float mesh_IoS(
    const std::vector<Point> &pts1, 
    const std::vector<Point> &pts2,
    const std::vector<std::vector<size_t>> &triangles1, 
    const std::vector<std::vector<size_t>> &triangles2,
    const float smallest_interval = 1.0,
    const int num_ray_per_pt = 2)
{
    // 0. create triangle set
    // std::list<Triangle> triangle_set1, triangle_set2;
    // for (const auto &t : triangles1)
    // {
    //     Triangle tri(pts1[t[0]], pts1[t[1]], pts1[t[2]]);
    //     if (tri.is_degenerate()) continue;
    //     triangle_set1.push_back(Triangle(pts1[t[0]], pts1[t[1]], pts1[t[2]]));
    // }
    // for (const auto &t : triangles2)
    // {
    //     Triangle tri(pts2[t[0]], pts2[t[1]], pts2[t[2]]);
    //     if (tri.is_degenerate()) continue;
    //     triangle_set2.push_back(Triangle(pts2[t[0]], pts2[t[1]], pts2[t[2]]));
    // }
        
    // 1. Calculate the bounds of two meshes and create two grids
    BBox bbox1 = CGAL::bounding_box(pts1.begin(), pts1.end()),
         bbox2 = CGAL::bounding_box(pts2.begin(), pts2.end());
    
    // if bbox not overlap, then return 0 directly
    if (bbox1.xmin() >= bbox2.xmax() || bbox2.xmin() >= bbox1.xmax() ||
        bbox1.ymin() >= bbox2.ymax() || bbox2.ymin() >= bbox1.ymax() ||
        bbox1.zmin() >= bbox2.zmax() || bbox2.zmin() >= bbox1.zmax())
        return 0;

    float xmin = std::min(CGAL::to_double(bbox1.xmin()), CGAL::to_double(bbox2.xmin())),
          xmax = std::max(CGAL::to_double(bbox1.xmax()), CGAL::to_double(bbox2.xmax())),
          ymin = std::min(CGAL::to_double(bbox1.ymin()), CGAL::to_double(bbox2.ymin())),
          ymax = std::max(CGAL::to_double(bbox1.ymax()), CGAL::to_double(bbox2.ymax())),
          zmin = std::min(CGAL::to_double(bbox1.zmin()), CGAL::to_double(bbox2.zmin())),
          zmax = std::max(CGAL::to_double(bbox1.zmax()), CGAL::to_double(bbox2.zmax()));
    
    float diag = std::sqrt((xmax - xmin) * (xmax - xmin) + 
                           (ymax - ymin) * (ymax - ymin) + 
                           (zmax - zmin) * (zmax - zmin));
    float diag1 = std::sqrt(CGAL::to_double(
                            (bbox1.xmax() - bbox1.xmin()) * (bbox1.xmax() - bbox1.xmin())
                          + (bbox1.ymax() - bbox1.ymin()) * (bbox1.ymax() - bbox1.ymin())
                          + (bbox1.zmax() - bbox1.zmin()) * (bbox1.zmax() - bbox1.zmin())));
    float diag2 = std::sqrt(CGAL::to_double(
                            (bbox2.xmax() - bbox2.xmin()) * (bbox2.xmax() - bbox2.xmin())
                          + (bbox2.ymax() - bbox2.ymin()) * (bbox2.ymax() - bbox2.ymin())
                          + (bbox2.zmax() - bbox2.zmin()) * (bbox2.zmax() - bbox2.zmin())));
    float min_diag = std::min(diag1, diag2);
    float interval  = min_diag / 10;
    interval = std::max(interval, smallest_interval);
    
    int X = int(std::ceil(xmax - xmin) / interval),
        Y = int(std::ceil(ymax - ymin) / interval),
        Z = int(std::ceil(zmax - zmin) / interval);
    size_t num_sampled = X * Y * Z;
    
    SDF_Points vtx1(pts1.size(), 3), vtx2(pts2.size(), 3), sample_pts(num_sampled, 3);
    SDF_Triangles tri1(triangles1.size(), 3), tri2(triangles2.size(), 3);
    for (int i = 0; i < pts1.size(); i++)
    {
        Eigen::RowVector3f pt {{static_cast<float>(CGAL::to_double(pts1[i].x())), 
                                static_cast<float>(CGAL::to_double(pts1[i].y())), 
                                static_cast<float>(CGAL::to_double(pts1[i].z()))}};
        vtx1.row(i).noalias() = pt;
    }
    for (int i = 0; i < pts2.size(); i++)
    {
        Eigen::RowVector3f pt {{static_cast<float>(CGAL::to_double(pts2[i].x())), 
                                static_cast<float>(CGAL::to_double(pts2[i].y())), 
                                static_cast<float>(CGAL::to_double(pts2[i].z()))}};
        vtx2.row(i).noalias() = pt;
    }

    int idx = 0;
    for (const auto &t : triangles1)
    {    
        Eigen::Matrix<std::uint32_t, 1, 3> index {{
            static_cast<std::uint32_t>(t[0]), static_cast<std::uint32_t>(t[1]), static_cast<std::uint32_t>(t[2])}}; 
        tri1.row(idx).noalias() = index; 
        idx++; 
    }
    idx = 0;
    for (const auto &t : triangles2)
    {       
        Eigen::Matrix<std::uint32_t, 1, 3> index {{
            static_cast<std::uint32_t>(t[0]), static_cast<std::uint32_t>(t[1]), static_cast<std::uint32_t>(t[2])}}; 
        tri2.row(idx).noalias() = index; 
        idx++; 
    }

    for (int i = 0; i < X; i++)
    {
        for (int j = 0; j < Y; j++)
        {
            for (int k = 0; k < Z; k++)
            {
                Eigen::RowVector3f pt {{xmin + i * interval, ymin + j * interval, zmin + k * interval}};
                sample_pts.row(i * Y * Z + j * Z + k).noalias() = pt;
            }
        }
    }

    // #if DEBUG_SYMMETRY == 1
    // std::cout << "xmin=" << xmin << ", xmax=" << xmax << "\n"
    //           << "ymin=" << ymin << ", ymax=" << ymax << "\n"
    //           << "zmin=" << zmin << ", zmax=" << zmax << "\n";
    // std::cout << "diag=" << diag << ", diag1=" << diag1 << ", diag2=" << diag2 
    //           << ", min_diag=" << min_diag << ", interval: "<< interval << "\n"; 
    // std::cout << "#pts1=" << pts1.size() << ", #pts2=" << pts2.size() << "\n";
    // std::cout << "#tri1=" << triangles1.size() << ", #tri2=" << triangles2.size() << "\n";
    // std::cout << "X=" << X << ", Y=" << Y << ", Z=" << Z << "\n";
    // #endif
    
    sdf::SDF sdf1(vtx1, tri1), sdf2(vtx2, tri2);
    SDF_Result query1(num_sampled), query2(num_sampled);
    sdf1.contains(sample_pts, query1);
    sdf2.contains(sample_pts, query2);
    int volume_self = (query1.array() > 0).count(), volume_other = (query2.array() > 0).count();
    int volume_inter = ((query1.array() > 0) * (query2.array() > 0)).count();

 
    float ios = 0;
    if (volume_self == 0 || volume_other == 0) 
        ios = 0;
    else
    {
        ios = float(volume_inter) / volume_self;
        ios = std::clamp(ios, float(0.0), float(1.0));
    }

    #if DEBUG_SYMMETRY == 1

    std::cout << "mesh ios=" << ios 
              << " intersection=" << volume_inter 
              << " self=" << volume_self << "\n";
    #endif

    return ios;
}

Eigen::Vector3f get_centroid(
    const std::vector<Eigen::Vector3f> &vtxs)
{
    Eigen::Vector3f c = Eigen::Vector3f::Zero();
    for (const auto &v : vtxs)
        c += v;
    return c / vtxs.size();
}

void sample_points_in_cell3(
    const LCC_3 &lcc, Dart_const_descriptor &cell_dart, 
    const float unit_volume, const int min_pts, const int max_pts, 
    SDF_Points &query_points)
{
    float volume = volume_of_cell3(lcc, cell_dart);
    int num_pts = std::max(min_pts, std::min(max_pts, int(std::ceil(volume/unit_volume))));

    std::vector<Eigen::Vector3f> vtxs;
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::const_iterator
        v = lcc.one_dart_per_incident_cell<0, 3>(cell_dart).begin();
        v != lcc.one_dart_per_incident_cell<0, 3>(cell_dart).end(); ++v)
        vtxs.push_back(Eigen::Vector3f(
            static_cast<float>(CGAL::to_double(lcc.point(v).x())),
            static_cast<float>(CGAL::to_double(lcc.point(v).y())),
            static_cast<float>(CGAL::to_double(lcc.point(v).z()))));
    size_t num_vtx = vtxs.size();
    std::cout << "num_points: " << num_pts << ", num_vtx=" << num_vtx << '\n';
    Eigen::Vector3f centroid = get_centroid(vtxs);
    std::vector<float> dist_to_centroid(num_vtx);
    float max_dist = 0;
    for (size_t i = 0; i < num_vtx; i++)
    {
        dist_to_centroid[i] = (vtxs[i] - centroid).norm();
        if (dist_to_centroid[i] > max_dist)
            max_dist = dist_to_centroid[i];
    }
    std::vector<float> vtx_weights(num_vtx);
    for (size_t i = 0; i < num_vtx; i++)
        vtx_weights[i] = dist_to_centroid[i] / max_dist;
    
    Eigen::MatrixXf weights = Eigen::MatrixXf::Zero(num_pts + num_vtx, num_vtx);
    // Eigen::MatrixXf weights = random;
    // weights = weights.array() + 1.0f;
    // weights = weights.array() / (weights.rowwise().sum().replicate(1, num_vtx)).array();
    for (size_t i = 0; i < num_vtx; i++)
    {
        for (size_t j = 0; j < num_pts/2; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());  //here you could also set a seed
            float hi = (1 - weights.row(j).sum()) * vtx_weights[(j % num_vtx + i) % num_vtx];
            std::uniform_real_distribution<float> dis(0, hi);
            weights(j, (j % num_vtx + i) % num_vtx) = dis(gen);
        }
    }

    for (size_t i = 0; i < num_vtx; i++)
    {
        for (size_t j = num_pts/2; j < num_pts; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());  //here you could also set a seed
            std::uniform_real_distribution<float> dis(0, 1);
            weights(j, i) = dis(gen);
        }
    }

    query_points = SDF_Points(num_pts + num_vtx, 3);
    for (size_t i = 0; i < num_pts; i++)
    {
        weights.row(i) = weights.row(i).array() / weights.row(i).sum();
        Eigen::RowVector3f pt = Eigen::RowVector3f::Zero();
        for (size_t j = 0; j < num_vtx; j++)
            pt = pt.array() + vtxs[j].transpose().array() * weights(i, j);
        query_points.row(i) = pt;
    }
    float other_weight = 0.1 / (num_vtx-1);
    for (size_t i = 0; i < num_vtx; i++)
    {
        Eigen::RowVector3f pt = Eigen::RowVector3f::Zero();
        for (size_t j = 0; j < num_vtx; j++)
            if (j != i)
                pt = pt.array() + vtxs[j].transpose().array() * other_weight;
            else
                pt = pt.array() + vtxs[j].transpose().array() * 0.9;
        query_points.row(num_pts + i) = pt;
    }
}

float mesh_IoS(
    const sdf::SDF &other_sdf, const SDF_Points & query_points)
{
    SDF_Result query_result(query_points.rows());
    other_sdf.contains(query_points, query_result);
    int num_pts = query_points.rows();
    int num_in = query_result.array().count();
    float ios = float(num_in) / num_pts;
    return ios;
}

float IoS(
    const LCC_3 &lcc,
    Dart_const_descriptor cell_dart,
    const float sample_unit_volume,
    const int min_point, const int max_point)
{
    float ios = 0;
    // sample the query points
    SDF_Points query_points;
    sample_points_in_cell3(lcc, cell_dart, sample_unit_volume, min_point, max_point, query_points);
    std::vector<Point> pts_this;
    for (int i = 0; i < query_points.rows(); i++)
    {
        Eigen::RowVector3f pt = query_points.row(i);
        pts_this.push_back(Point(pt[0], pt[1], pt[2]));
    }
    BBox bbox_this = CGAL::bounding_box(pts_this.begin(), pts_this.end());

    for (LCC_3::One_dart_per_cell_range<3>::const_iterator
        it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        Volume_idx vid = lcc.info<3>(it).first;
        std::vector<Point> pts2;
        std::vector<std::vector<size_t>> triangles2;
        extract_triangle_mesh(lcc, vid, pts2, triangles2);
        BBox bbox2 = CGAL::bounding_box(pts2.begin(), pts2.end());
        if (bbox_this.xmin() >= bbox2.xmax() || bbox2.xmin() >= bbox_this.xmax() ||
            bbox_this.ymin() >= bbox2.ymax() || bbox2.ymin() >= bbox_this.ymax() ||
            bbox_this.zmin() >= bbox2.zmax() || bbox2.zmin() >= bbox_this.zmax())
            continue;

        SDF_Points vtx2(pts2.size(), 3);
        SDF_Triangles tri2(triangles2.size(), 3);
        for (int i = 0; i < pts2.size(); i++)
        {
            Eigen::RowVector3f pt {{static_cast<float>(CGAL::to_double(pts2[i].x())), 
                                    static_cast<float>(CGAL::to_double(pts2[i].y())), 
                                    static_cast<float>(CGAL::to_double(pts2[i].z()))}};
            vtx2.row(i) = pt;
        }
        for (int i = 0; i < triangles2.size(); i++)
        {    
            Eigen::Matrix<std::uint32_t, 1, 3> index {{
                static_cast<std::uint32_t>(triangles2[i][0]), 
                static_cast<std::uint32_t>(triangles2[i][1]), 
                static_cast<std::uint32_t>(triangles2[i][2])}}; 
            tri2.row(i) = index; 
        }
        sdf::SDF other_sdf(vtx2, tri2, true);
        float this_ios = mesh_IoS(other_sdf, query_points);
        if (this_ios > ios) ios = this_ios;
    }
    return ios;
}

bool is_cell3_on_one_side(
    const LCC_3 &lcc,
    Dart_const_descriptor cell_dart,
    const Plane* plane, 
    const Surface_mesh &mesh,
    const std::vector<Face_index> &faces,
    const float dist_thresh = 0.5,
    const float buff_dist = 0.5)
{
    bool on_pos = false, on_neg = false;
    std::vector<Point> pts_cell;
    // loop all vertices
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::const_iterator
        it=lcc.one_dart_per_incident_cell<0, 3>(cell_dart).begin(),
        itend=lcc.one_dart_per_incident_cell<0, 3>(cell_dart).end(); it!=itend; ++it)
    {
        pts_cell.push_back(lcc.point(it));
        if (plane->oriented_side(lcc.point(it)) == CGAL::ON_POSITIVE_SIDE)
        {
            float dist = dist_pt_to_plane(lcc.point(it), plane);
            if (dist > dist_thresh) on_pos = true;
        }
        else if (plane->oriented_side(lcc.point(it)) == CGAL::ON_NEGATIVE_SIDE)
        {
            float dist = dist_pt_to_plane(lcc.point(it), plane);
            if (dist > dist_thresh) on_neg = true;
        }
            
    }
    if (!on_pos && on_neg) return true;
    if (on_pos != on_neg)
        return true;
    else
    {
        std::vector<Point> pts_plane;
        for (const auto & f: faces)
        {
            Halfedge_index he = mesh.halfedge(f);
            do
            {
                Point pt = mesh.point(mesh.source(he));
                pts_plane.push_back(pt);
                he = mesh.next(he);
            } while (he != mesh.halfedge(f));
        }
        
        auto bbox_cell = CGAL::bounding_box(pts_cell.begin(), pts_cell.end());
        auto bbox_plane = CGAL::bounding_box(pts_plane.begin(), pts_plane.end());
        if (bbox_cell.xmin() >= bbox_plane.xmax() + buff_dist 
         || bbox_plane.xmin() >= bbox_cell.xmax() + buff_dist
         || bbox_cell.ymin() >= bbox_plane.ymax() + buff_dist 
         || bbox_plane.ymin() >= bbox_cell.ymax() + buff_dist
         || bbox_cell.zmin() >= bbox_plane.zmax() + buff_dist 
         || bbox_plane.zmin() >= bbox_cell.zmax() + buff_dist)
            return true;
        else
            return false;
    }
        
}

std::mt19937& get_rng() {
    // Safer seeding with time (random_device can be not availble)
    thread_local std::mt19937 rg{
        std::random_device{}() ^
        static_cast<unsigned int>(std::chrono::high_resolution_clock::now()
                                      .time_since_epoch()
                                      .count())};
    return rg;
}

bool are_plane_parallel(
    const Plane* plane1, const Plane* plane2, const float angle_thresh = 10)
{
    Eigen::Vector3d n1 (CGAL::to_double(plane1->a()), CGAL::to_double(plane1->b()), CGAL::to_double(plane1->c())),
                    n2 (CGAL::to_double(plane2->a()), CGAL::to_double(plane2->b()), CGAL::to_double(plane2->c()));
    n1.normalize();
    n2.normalize();
    float radian = std::acos(n1.dot(n2));
    float degree = radian * 180 / M_PI;
    if (degree < angle_thresh || degree > 180 - angle_thresh)
        return true;
    else
        return false;
}

size_t merge_small_cell(
    LCC_3 &lcc, const float& total_volume, const float min_vol_ratio)
{
    size_t merged_count = 0;
    bool merged = false; 
    float mean_vol = total_volume / lcc.number_of_attributes<3>();
    float small_thresh = total_volume * min_vol_ratio;
    std::cout << "merge small cells...\n";
    
    do{
        merged = false;
        std::cout << "new round\n";
        // mean_vol = total_volume / lcc.number_of_attributes<3>(); 
        // small_thresh = mean_vol * min_vol_ratio;
        Volume_idx this_v = -1, oppo_v = -1;
        // DO NOT MODIFY CM IN THE I-CELL LOOP 
        for (LCC_3::One_dart_per_cell_range<3>::iterator
            it=lcc.one_dart_per_cell<3>().begin(),
            itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
        {
            // Your logic here
            float volume = volume_of_cell3(lcc, it);
            if (volume < small_thresh)
            {
                this_v = lcc.info<3>(it).first;
                oppo_v = this_v;
                std::map<Volume_idx, float> adj_areas;
                for (LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator
                    f=lcc.one_dart_per_incident_cell<2, 3>(it).begin(),
                    fend=lcc.one_dart_per_incident_cell<2, 3>(it).end(); f!=fend; ++f)
                {
                    std::cout << "find face\n";
                    if (lcc.is_free<3>(f)) continue;
                    Dart_const_descriptor oppo_f = lcc.beta<3>(f);
                    oppo_v = lcc.info<3>(oppo_f).first;
                    bool oppo_in = lcc.info<3>(oppo_f).second;
                    
                    if (this_v != oppo_v && oppo_in) 
                    {
                        if (adj_areas.find(oppo_v) != adj_areas.end())
                            adj_areas[oppo_v] += area_of_cell2(lcc, f);
                        else
                            adj_areas[oppo_v] = area_of_cell2(lcc, f);
                    }
                }
                if (adj_areas.size() != 0) 
                {
                    merged = true;
                    // find the adjacent cell with largest area to merge
                    float largest = 0;
                    for (auto [vid, area] : adj_areas)
                    {   
                        if (area > largest)
                        {
                            largest = area;
                            oppo_v = vid;
                        }
                    }
                    break;
                }
            }
            
        }
        if (merged)
        {
            merged_count++;
            std::vector<Dart_descriptor> edges_to_remove;
            Dart_descriptor this_cell_dart;
            for (LCC_3::One_dart_per_cell_range<3>::iterator
                vit=lcc.one_dart_per_cell<3>().begin(),
                vitend=lcc.one_dart_per_cell<3>().end(); vit!=vitend; ++vit)
            {
                if (lcc.info<3>(vit).first == this_v)
                {
                    this_cell_dart = vit;
                    break;
                }
            }
            for (LCC_3::One_dart_per_incident_cell_range<1, 3>::iterator
                    e=lcc.one_dart_per_incident_cell<1, 3>(this_cell_dart).begin(),
                    eend=lcc.one_dart_per_incident_cell<1, 3>(this_cell_dart).end(); e!=eend; ++e)
            {
                if (lcc.is_free<3>(e) || lcc.is_free<2>(e)) continue;
                if (lcc.info<3>(lcc.beta<3>(e)).first == oppo_v)
                {
                    std::set<const Plane*> planes;
                    for (LCC_3::One_dart_per_incident_cell_range<2, 1>::iterator
                        f=lcc.one_dart_per_incident_cell<2, 1>(e).begin(),
                        fend=lcc.one_dart_per_incident_cell<2, 1>(e).end(); f!=fend; ++f)
                        planes.insert(lcc.info<2>(f).second);
                    if (planes.size() == 1)
                        edges_to_remove.push_back(e);
                }
            }
            std::cout << "collect edges=" << edges_to_remove.size() << "\n";
            for (auto e : edges_to_remove)
            {
                if (lcc.is_removable<1>(e))
                {
                    lcc.remove_cell<1>(e);
                    std::cout << "removed edge\n";
                }
                else
                    std::cout << "unremovable edge\n";
            }

            std::vector<Dart_descriptor> faces_to_remove;
            for (LCC_3::One_dart_per_cell_range<3>::iterator
                vit=lcc.one_dart_per_cell<3>().begin(),
                vitend=lcc.one_dart_per_cell<3>().end(); vit!=vitend; ++vit)
            {
                if (lcc.info<3>(vit).first == this_v)
                {
                    this_cell_dart = vit;
                    break;
                }
            }
            for (LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator
                    ff=lcc.one_dart_per_incident_cell<2, 3>(this_cell_dart).begin(),
                    ffend=lcc.one_dart_per_incident_cell<2, 3>(this_cell_dart).end(); ff!=ffend; ++ff)
            {
                if (lcc.is_free<3>(ff) || lcc.is_free<2>(ff)) continue;
                if (lcc.info<3>(lcc.beta<3>(ff)).first == oppo_v)
                    faces_to_remove.push_back(ff);
            }
            std::cout << "collect faces=" << faces_to_remove.size() << "\n";
            for (auto ff : faces_to_remove)
            {
                if (lcc.is_removable<2>(ff))
                {
                    lcc.remove_cell<2>(ff);
                    std::cout << "removed face\n";
                }
                else
                    std::cout << "unremovable face\n";
                
            }
            std::cout << "update index of cm ...\n";
            Vertex_idx vtx_id = 0;
            for(LCC_3::One_dart_per_cell_range<0>::iterator
                it= lcc.one_dart_per_cell<0>().begin(),
                itend= lcc.one_dart_per_cell<0>().end(); it!=itend; ++it)
                lcc.info<0>(it) = vtx_id++;
            
            Edge_idx eid = 0;
            for(LCC_3::One_dart_per_cell_range<1>::iterator
                it= lcc.one_dart_per_cell<1>().begin(),
                itend= lcc.one_dart_per_cell<1>().end(); it!=itend; ++it)
                lcc.info<1>(it) = eid++;

            Face_idx fid = 0;
            for(LCC_3::One_dart_per_cell_range<2>::iterator
                it= lcc.one_dart_per_cell<2>().begin(),
                itend= lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
            {
                lcc.info<2>(it).first = fid++;
            }

            Volume_idx vid = 0;
            for(LCC_3::One_dart_per_cell_range<3>::iterator
                it= lcc.one_dart_per_cell<3>().begin(),
                itend= lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
                    lcc.info<3>(it).first = vid++;
            std::cout << "updated!\n";
        }
    }while(merged);

    return merged_count;
}