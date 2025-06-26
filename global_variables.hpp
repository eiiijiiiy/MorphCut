#pragma once

#include "util.hpp"

class Component;
struct Evaluate_weights;

namespace Global_var
{
    sdf::SDF sdf;
    bool use_raycasting;
    Tree tree;
    std::unordered_map<Dart_const_descriptor, std::pair<SDF_Points, SDF_Result>> cell_sdf;
    std::string input_file;
    float area_horizontal = 0;
    float concave_thresh;
    float cell_size;
    size_t min_points;
    size_t max_points;
    float sampled_unit_volume;
    float surface_area;
    float min_cut_factor;
    float same_vol_factor;
    float mv_weight;
    float sym_weight;
    float cen_weight;
    float de_weight;
    int max_children;
    Eigen::Array<bool, Eigen::Dynamic, 1> vote_map;
    Surface_mesh input_mesh;
    bool use_repair;
    Surface_mesh repaired_mesh;
    BBox box;
    LCC_3 lcc;
    std::unordered_map<int, const Plane*> face_planes; 
    std::unordered_map<Dart_descriptor, In_cell_pts> in_cell_pts_map;
    std::unordered_map<int, std::pair<bool, Dart_descriptor>> is_in_cell;
    std::vector<Point> grid_points;
    std::vector<const Plane*> rg_planes;
    std::unordered_map<const Plane*, Eigen::ArrayXi> on_rg_plane_pos;
    std::unordered_map<const Plane*, std::vector<Face_index>> plane_regions;
    std::unordered_map<const Plane*, const std::vector<Point>> plane_extended_convex_hulls;
    std::unordered_map<Dart_descriptor, std::vector<std::pair<const Plane*, float>>> to_cut_planes;
    int total_in;
    float angle_thresh;
    float radian_thresh;
    float dist_thresh;
    // float mcts_c = 0.01;
    float total_volume;
    std::map<int, Dart_descriptor> init_comp;
    std::unordered_map<const Plane*, Eigen::ArrayXi> on_imp_plane_pos_cell;
    std::set<const Plane*> all_face_planes;
    std::set<const Plane*> all_face_planes_copy;
    std::vector<float> volumes_of_int_cells;
    float max_volume_deviation;
    Evaluate_weights eval_weights;
    int num_init_cells;
}