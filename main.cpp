#include "util.hpp"
#include "global_variables.hpp"
#include "reconstruction.hpp"
#include "Node.hpp"
#include "Component.hpp"
#include "specific_search.hpp"
// #include "bnb.hpp"
float Node::cvx_thresh_;
std::function<float(const Node*)> Node::m_g;
std::function<float(const Node*)> Node::m_f;
float Component::cvx_thresh_;
int Component::max_try_ = 0;
enum Cut_mode Component::m_s_cut_mode;
bool Component::m_with_bisector;
float Component::m_min_cut_factor;

int main(int argc, char * argv[])
{
    #if LOG_TIME == 1
    auto start_prep = std::chrono::steady_clock::now();
    #endif
    Global_var::input_file = std::string(argv[1]);
    
    std::cout << "input_file: " << Global_var::input_file << "\n";

    /***************** Read & translate mesh *******************/
    if (!PMP::IO::read_polygon_mesh(Global_var::input_file, Global_var::input_mesh)) {
        std::cerr << "ERROR: cannot read the input file:"
            << Global_var::input_file << '\n';
        return EXIT_FAILURE;
    }
    std::string repaired_file = Global_var::input_file + "_repaired.obj";
    Global_var::use_repair = true;
    if (!PMP::IO::read_polygon_mesh(repaired_file, Global_var::repaired_mesh)) 
        Global_var::use_repair = false;
    
    std::cout << "use repair: " << Global_var::use_repair << "\n";

    Vector original_center;
    original_center = move_mesh_to_origin(Global_var::input_mesh);
    if (Global_var::use_repair)
        move_mesh(Global_var::repaired_mesh, original_center);
    
    std::cout << "original center: " << original_center << '\n';
    if (!CGAL::is_triangle_mesh(Global_var::input_mesh)) {
        std::cout << "input geometry is not triangulated. Triangulating ...\n";
        PMP::triangulate_faces(Global_var::input_mesh);
        std::cout << "done\n";
    }
    std::cout << "input mesh now is triangulated\n";
    Global_var::box = CGAL::bounding_box(
                Global_var::input_mesh.points().begin(), Global_var::input_mesh.points().end());
    float X = CGAL::to_double(Global_var::box.xmax() - Global_var::box.xmin()),
          Y = CGAL::to_double(Global_var::box.ymax() - Global_var::box.ymin()),
          Z = CGAL::to_double(Global_var::box.zmax() - Global_var::box.zmin());

    rj::Document config_doc;
    std::string config_file = Global_var::input_file + "_config.json";
    read_config(config_file.c_str(), config_doc);
    std::cout << "loaded config file: " << config_file << '\n';
    // float cell_size = std::sqrt(X*X + Y*Y + Z*Z) * 0.05 * factor_cell_size;
    // double max_distance = cell_size * 0.05;
    // Global_var::concave_thresh = 1.2;
    // double max_angle = 10.0;
    // Global_var::angle_thresh = std::max(5.0, max_angle/5.0);
    // Global_var::dist_thresh = std::max(0.5, max_distance/5.0);
    // LCC::vol_thresh = CGAL::to_double(X * Y * Z) * 1e-3;
    float max_distance = config_doc["max_distance"].GetFloat();
    float max_angle = config_doc["max_angle"].GetFloat();
    size_t min_region = config_doc["min_region"].GetInt();
    Global_var::angle_thresh = config_doc["angle_thresh"].GetFloat();
    Global_var::radian_thresh = Global_var::angle_thresh / 180 * M_PI;
    Global_var::dist_thresh = config_doc["dist_thresh"].GetFloat();
    float cell_size = config_doc["cell_size"].GetFloat();
    Global_var::cell_size = cell_size;
    Global_var::min_points = config_doc["min_points"].GetInt();
    Global_var::max_points = config_doc["max_points"].GetInt();
    float sigma = config_doc["sigma"].GetFloat();
    float cvx_thresh = 0.999;
    if (config_doc.HasMember("cvx_thresh"))
        cvx_thresh = config_doc["cvx_thresh"].GetFloat();
    float extended_plane_ch;
    if (config_doc.HasMember("extended_plane_ch"))
        extended_plane_ch = config_doc["extended_plane_ch"].GetFloat();
    else
        extended_plane_ch = 3 * cell_size;
    Global_var::concave_thresh = config_doc["concave_thresh"].GetFloat();
    // if (config_doc.HasMember("min_cut_factor"))
    //     Global_var::min_cut_factor = config_doc["min_cut_factor"].GetFloat();
    // else
    //     Global_var::min_cut_factor = 0.001;
    if (config_doc.HasMember("same_vol_factor"))
        Global_var::same_vol_factor = config_doc["same_vol_factor"].GetFloat();
    else
        Global_var::same_vol_factor = 0.001;
    if (config_doc.HasMember("max_children"))
        Global_var::max_children = config_doc["max_children"].GetInt();
    else
        Global_var::max_children = -1;

    std::cout << "configuration: "
              << " \t initial region: " << min_region << ", max_distance: " << max_distance << ", max_angle: " << max_angle << '\n'
              << " \t refined region: " << Global_var::angle_thresh << ", " << Global_var::dist_thresh << '\n'
              << " \t cell_size: " << cell_size << ", X: " << X << ", Y: " << Y << ", Z: " << Z << '\n'
              << " \t extended_plane_ch: " << extended_plane_ch << '\n'
              << " \t refined region: " << Global_var::angle_thresh << ", " << Global_var::dist_thresh << '\n'
              << " \t concave_thresh: " << Global_var::concave_thresh << '\n';
            //   << " \t MCTS: max_iter: " << max_iter << ", max_depth: " << max_depth << '\n';


    std::cout << "calculating areas of faces\n";
    calculate_areas_faces_3d(Global_var::input_mesh);
    std::optional<Surface_mesh::Property_map<Face_index, float>> face_areas = 
        Global_var::input_mesh.property_map<Face_index, float>("f:area");
    assert(face_areas.has_value());
    /***************** Region growing *******************/
    Neighbor_query query(Global_var::input_mesh);
    Region_type region_type(
        Global_var::input_mesh, CGAL::parameters::
        maximum_distance(FT(max_distance)).
        maximum_angle(FT(max_angle)).
        minimum_region_size(min_region));

    Sorting sorting(Global_var::input_mesh, query);
    sorting.sort();

    const auto& face_range = faces(Global_var::input_mesh);
    std::cout << "* number of input faces: " << face_range.size() << '\n';

    Region_growing region_growing(
        face_range, sorting.ordered(), query, region_type);

    std::vector<Primitive_and_region> regions;
    region_growing.detect(back_inserter(regions));
    std::cout << "* number of detected regions: " << regions.size() << '\n';

    // CGAL::IO::write_PLY(rg_output, Global_var::input_mesh);
    /***************** Pre-locate grid points of planes *******************/
    std::unordered_map<const Plane*, float> plane_area_ratio;
    std::vector<std::pair<const Plane*, float>> to_cut_planes_of_first_cell;

    for (const auto &region : regions)
    {
        const Plane * normalized_p = normalize_plane(&region.first);
        auto [dup, ext_plane] = is_duplicate_plane(Global_var::input_mesh,
            normalized_p, region.second, Global_var::plane_regions, 
            Global_var::angle_thresh, Global_var::dist_thresh); // 有bug？
        // dup = false; 
        if (!dup)
        {
            Global_var::rg_planes.push_back(normalized_p);
            for (auto f : region.second)
                Global_var::plane_regions[normalized_p].push_back(f);
        }
        else
        {
            for (auto f : region.second)
                Global_var::plane_regions[ext_plane].push_back(f);
        }
    }
    // filter out small planar regions
    std::vector<const Plane*> to_remove;
    // FT cell_size (std::stod(std::string(argv[5])));

    for (auto plane : Global_var::rg_planes)
    {
        float area = 0;
        for (auto f : Global_var::plane_regions[plane])
            area += face_areas.value()[f];
        
        if (area < cell_size * cell_size) 
            to_remove.push_back(plane);
        else
            to_cut_planes_of_first_cell.push_back(std::make_pair(plane, area));
    }
    for (auto plane : to_remove)
    {
        std::vector<const Plane*>::iterator it = std::find(Global_var::rg_planes.begin(), Global_var::rg_planes.end(), plane);
        Global_var::rg_planes.erase(it);
        Global_var::plane_regions.erase(plane);
    }
    std::cout << "* number of region growing planes (refined): " << Global_var::rg_planes.size() << '\n';
    bool created;
    Surface_mesh::Property_map<Face_index, CGAL::IO::Color> f_colors;
    boost::tie(f_colors, created) =
        Global_var::input_mesh.add_property_map<Face_index, CGAL::IO::Color>(
            "f:color", CGAL::IO::Color(CGAL::IO::gray()));
    
    int color_idx = 0;
    for (auto plane : Global_var::rg_planes)
    {
        auto rg_color = pastelColorMap[color_idx % pastelColorMap.size()];
        // std::cout << "plane: " << plane << " color: " << rg_color[0] << ' ' << rg_color[1] << ' ' << rg_color[2] << '\n';
        for (auto f : Global_var::plane_regions[plane])
            f_colors[f] = CGAL::IO::Color(rg_color[0], rg_color[1], rg_color[2]);
        color_idx++;
    }
    

    compute_extended_convex_hull_of_planes(
            Global_var::plane_regions, Global_var::plane_extended_convex_hulls, 
            Global_var::input_mesh, extended_plane_ch);

    /***************** prepare for SDF *******************/
    std::map<Vertex_index, size_t> vidx_map;
    size_t idx = 0;
    Surface_mesh mesh;
    if (Global_var::use_repair)
        mesh = Global_var::repaired_mesh;
    else
        mesh = Global_var::input_mesh;

    size_t num_pts = num_vertices(mesh);
    SDF_Points pts(num_pts, 3);
    
    for (Vertex_index v : vertices(mesh))
    {
        Point p = mesh.point(v);
        Eigen::RowVector3f pt {{
                        static_cast<float>(CGAL::to_double(p.x())), 
                        static_cast<float>(CGAL::to_double(p.y())), 
                        static_cast<float>(CGAL::to_double(p.z()))}};
        pts.row(idx) = pt;
        vidx_map[v] = idx;
        idx++;
    }

    size_t num_faces = std::distance(
                faces(mesh).first, 
                faces(mesh).second);
    SDF_Triangles tris(num_faces, 3);
    idx = 0;
    for (Face_index f: faces(mesh))
    {
        Halfedge_index e = mesh.halfedge(f);
        Eigen::Matrix<std::uint32_t, 1, 3> index {{
            static_cast<std::uint32_t>(vidx_map[mesh.source(e)]),
            static_cast<std::uint32_t>(vidx_map[mesh.source(mesh.next(e))]),
            static_cast<std::uint32_t>(vidx_map[mesh.source(mesh.next(mesh.next(e)))])
        }};
        tris.row(idx) = index;
        idx++;
    }
    Global_var::sdf.init(pts, tris, true);
    Global_var::use_raycasting = false;
    std::cout << "use repair: " << Global_var::use_repair << "\n";
    if (Global_var::use_repair)
        Global_var::tree = Tree(faces(Global_var::repaired_mesh).first, 
                    faces(Global_var::repaired_mesh).second, Global_var::repaired_mesh);
    else
        Global_var::tree = Tree(faces(Global_var::input_mesh).first, 
                    faces(Global_var::input_mesh).second, Global_var::input_mesh);
    
    Global_var::area_horizontal = CGAL::to_double(
        Global_var::box.xmax() - Global_var::box.xmin()) * CGAL::to_double(Global_var::box.ymax() - Global_var::box.ymin());
    
    Global_var::lcc.onsplit_functor<1>()=Split_functor_1(Global_var::lcc);
    Global_var::lcc.onsplit_functor<2>()=Split_functor_2(Global_var::lcc);
    Global_var::lcc.onsplit_functor<3>()=Split_functor_3(Global_var::lcc);
    Global_var::lcc.onmerge_functor<3>() = Merge_functor_3();
    
    Dart_descriptor first_cell = box_cell();
    // size_t num_grid_points = Global_var::grid_points.size();
    // std::cout << "* number of sampled points " << num_grid_points << '\n';

    // /***************** Vote in/out *******************/
    // const size_t num_ray_per_pt = 5;
    // const float in_thresh = 0.5;
    // vote_in_out(num_ray_per_pt, in_thresh);
    // std::cout << "* number of the in points: " << Global_var::vote_map.count() << '\n';
    // std::string in_pts_file = Global_var::input_file + "_in.xyz";
    // std::string out_pts_file = Global_var::input_file + "_out.xyz";
    // std::ofstream in_pts_o(in_pts_file);
    // std::ofstream out_pts_o(out_pts_file);
    // std::vector<Point> shirnk_pts;
    // for (size_t i = 0; i < num_grid_points; i++)
    // {
    //     if (Global_var::vote_map(i))
    //     {
    //         in_pts_o << Global_var::grid_points[i].x() << ' ' 
    //                  << Global_var::grid_points[i].y() << ' ' 
    //                  << Global_var::grid_points[i].z() << '\n';
    //         shirnk_pts.push_back(Global_var::grid_points[i]);
    //     }
    //     else
    //         out_pts_o << Global_var::grid_points[i].x() << ' ' 
    //                   << Global_var::grid_points[i].y() << ' ' 
    //                   << Global_var::grid_points[i].z() << '\n';
    // }
    // Global_var::grid_points = shirnk_pts;
    // num_grid_points = Global_var::grid_points.size();
    // in_pts_o.close();
    // out_pts_o.close();

    Global_var::to_cut_planes.insert({first_cell, to_cut_planes_of_first_cell});
    // pre_locate_grid_points_of_planes();

    // Global_var::in_cell_pts_map[first_cell] = std::vector<int>(num_grid_points);
    std::iota(Global_var::in_cell_pts_map[first_cell].begin(), 
              Global_var::in_cell_pts_map[first_cell].end(), 0);
    Global_var::sampled_unit_volume = cell_size * cell_size * cell_size;

    #if LOG_TIME == 1
    std::chrono::duration<float> duration_prep = std::chrono::steady_clock::now() - start_prep;
    std::cout << "[Duration] preparation: " << duration_prep.count() << " s\n";
    std::ofstream time_o(Global_var::input_file + "_time.txt");
    time_o << "preparation: " << duration_prep.count() << " s\n";
    time_o.close();
    #endif

    std::string rg_output = Global_var::input_file + "_rg.off";
    CGAL::IO::write_OFF(rg_output, Global_var::input_mesh, CGAL::parameters::face_color_map(f_colors));
    
    /**************** Iteratively cut *******************/
    #if LOG_TIME == 1
    auto start_reco = std::chrono::steady_clock::now();
    #endif
    
    iterative_cut_for_reconstruction();
    std::cout << "done reconstruction \n";
    
    #if LOG_TIME == 1
    std::chrono::duration<float> duration_reco = std::chrono::steady_clock::now() - start_reco;
    std::cout << "[Duration] reconstruction: " << duration_reco.count() << " s\n";
    time_o.open(Global_var::input_file + "_time.txt", std::ios::app);
    time_o << "reconstruction: " << duration_reco.count() << " s\n";
    time_o.close();
    #endif

    // /**************** Decompose *******************/
    std::string cut_mode = config_doc.HasMember("cut_mode") ? config_doc["cut_mode"].GetString() : "A";
    std::cout << "cut mode = " << cut_mode << "\n";
    Component::m_s_cut_mode = 
        cut_mode == "HV" ? Cut_mode::HV : 
        (cut_mode == "H" ? Cut_mode::H : 
        (cut_mode == "HI" ? Cut_mode::HI : 
        (cut_mode == "V" ? Cut_mode::V : 
        (cut_mode == "VI" ? Cut_mode::VI : Cut_mode::A))));
    
    if (config_doc.HasMember("with_bisector"))
        Component::m_with_bisector = config_doc["with_bisector"].GetBool();
    else
        Component::m_with_bisector = true;

    switch (Component::m_s_cut_mode)
    {
        case Cut_mode::HV:
            std::cout << "HV CUT\n";
            break;
        case Cut_mode::H:
            std::cout << "H CUT \n";
            break;
        default:
            std::cout << "other cut\n";
            break;
    }

    std::vector<float> lod_threshs;
    if (config_doc.HasMember("min_cut_factor"))
    {
        lod_threshs.clear();
        assert(config_doc["min_cut_factor"].IsArray());
        for (auto &v : config_doc["min_cut_factor"].GetArray())
            lod_threshs.push_back(v.GetFloat());
    }
    else
        lod_threshs = {0.1, 0.01, 0.001};
    

    
    if (config_doc.HasMember("max_vol_weight"))
        Global_var::mv_weight = config_doc["max_vol_weight"].GetFloat();
    else
        Global_var::mv_weight = 1.0;
    
    if (config_doc.HasMember("symmetry_weight"))
        Global_var::sym_weight = config_doc["symmetry_weight"].GetFloat();
    else
        Global_var::sym_weight = 1.0;
    
    if (config_doc.HasMember("central_weight"))
        Global_var::cen_weight = config_doc["central_weight"].GetFloat();
    else
        Global_var::cen_weight = 1.0;
    
    if (config_doc.HasMember("decentral_weight"))
        Global_var::de_weight = config_doc["decentral_weight"].GetFloat();
    else
        Global_var::de_weight = 1.0;
    
    std::string param_file = Global_var::input_file + "_param.txt";
    std::ofstream param_o(param_file);
    param_o << "origional_center:(" << original_center.x() << "," << original_center.y() << "," << original_center.z() << ")\n"
            << "planar_region=" << Global_var::rg_planes.size() << '\n'
            << "max_distance=" << max_distance << ", max_angle=" << max_angle << ", min_region=" << min_region << '\n'
            << "angle_thresh=" << Global_var::angle_thresh << ", dist_thresh=" << Global_var::dist_thresh << '\n'
            << "cell_size=" << cell_size << ", concave_thresh=" << Global_var::concave_thresh << '\n'
            << "min_points=" << Global_var::min_points << ", max_points=" << Global_var::max_points << '\n'
            << "sigma=" << sigma << ", max_children=" << Global_var::max_children << ", same_vol_factor=" << Global_var::same_vol_factor << '\n';
    param_o << "lod_threshs=";
    for (auto v : lod_threshs)
        param_o << v << ' ';
    param_o << '\n';

    bool fewer_comp = config_doc.HasMember("fewer_comp") && config_doc["fewer_comp"].GetBool(),
         max_vol = config_doc.HasMember("max_vol") && config_doc["max_vol"].GetBool(),
         symmetry = config_doc.HasMember("symmetry") && config_doc["symmetry"].GetBool(),
         central = config_doc.HasMember("central") && config_doc["central"].GetBool(),
         decentral = config_doc.HasMember("decentral") && config_doc["decentral"].GetBool();

    param_o << "Fewer_comp=" << fewer_comp << ", weight=" << 1.0 << '\n';
    param_o << "Max_vol=" << max_vol << ", weight=" << Global_var::mv_weight << '\n';
    param_o << "Symmetry=" << symmetry << ", weight=" << Global_var::sym_weight << '\n';
    param_o << "Central=" << central << ", weight=" << Global_var::cen_weight << '\n';
    param_o << "Decentral=" << decentral << ", weight=" << Global_var::de_weight << '\n';

    Component::m_min_cut_factor = 0.00001;
    Component* root_comp = new Component(0, Global_var::lcc, true);
    size_t num_cc_edges = root_comp->num_all_cc_edges();
    size_t num_cuts = root_comp->num_potential_cuts();
    param_o << "num_cc_edges=" << num_cc_edges << ", num_cuts=" << num_cuts << '\n';
    param_o.close();

    Node::m_f = f_general_n;
    std::function<float(const Node*)> g =
        [max_vol, symmetry, central, decentral](const Node* node) -> float
    {
        float value = float(node->num_comps()) / node->max_possible_comps()
           + node->reduced_delta_ch() / node->initial_delta_ch();
        if (max_vol)
            value += Global_var::mv_weight * (1 - node->max_volume() / node->sum_volume());
        if (symmetry)
            value += Global_var::sym_weight * node->convex_sym_iou();
        if (central)
            value += Global_var::cen_weight * node->get_ratio_comp_degree_2_stable();
        if (decentral)
            value += Global_var::de_weight * node->get_ratio_comp_more_than_degree_2_stable();
        return value;
    };
    Node::m_g = g;
    
    #if LOG_TIME == 1
    auto start_decomp = std::chrono::steady_clock::now();
    #endif

    Global_var::all_face_planes.clear();
    Global_var::all_face_planes = Global_var::all_face_planes_copy;
    std::string name = Global_var::input_file;
    if (fewer_comp) name += "_fewer_comp";
    if (max_vol) name += "_max_vol";
    if (symmetry) name += "_symmetry";
    if (central) name += "_central";
    if (decentral) name += "_decentral";
    
    if (lod_threshs.size() > 1)
    {
        std::vector<std::shared_ptr<const Node>> lod_results;
        if (symmetry)
        {
            std::vector<Point> pts;
            std::vector<std::vector<size_t>> triangles;
            extract_triangle_mesh(Global_var::lcc, 0, pts, triangles);
            std::pair<float, float> sym_axis = detect_symmetry(pts, triangles);
            bnb_lod(Global_var::lcc, lod_threshs, lod_results, sigma, name, sym_axis);
        }
        else
            bnb_lod(Global_var::lcc, lod_threshs, lod_results, sigma, name);
    }
    else // single lod decomposition, output all optimal results of the lod
    {
        std::vector<std::shared_ptr<const Node>> all_results;
        LCC_3 cur_lcc = Global_var::lcc;
        cur_lcc.onsplit_functor<1>() = Split_functor_1(cur_lcc);    
        cur_lcc.onsplit_functor<2>() = Split_functor_2(cur_lcc);
        cur_lcc.onsplit_functor<3>() = Split_functor_3(cur_lcc);
        cur_lcc.onmerge_functor<2>() = Merge_functor_2();
        cur_lcc.onmerge_functor<3>() = Merge_functor_3();
        Component* comp = new Component(0, cur_lcc, true);
        size_t num_cc_edges = comp->num_all_cc_edges();
        float volume_ch = comp->volume_ch();
        Node::cvx_thresh_ = 0.999;
        Component::cvx_thresh_ = 0.999;
        Component::m_min_cut_factor = lod_threshs[0];
        std::string search_name = name + "_search.txt";
        if (symmetry)
        {
            std::vector<Point> pts;
            std::vector<std::vector<size_t>> triangles;
            extract_triangle_mesh(Global_var::lcc, 0, pts, triangles);
            std::pair<float, float> sym_axis = detect_symmetry(pts, triangles);
            Node* root = new Node({comp}, cur_lcc, num_cc_edges + 1, 
                                    volume_ch, 0, sigma, sym_axis);
            
            bnb(root, all_results, search_name);
        }
        else
        {
            Node* root = new Node({comp}, cur_lcc, num_cc_edges + 1, volume_ch, 0, sigma);
            bnb(root, all_results, search_name);
        }

        for (size_t i = 0; i < all_results.size(); i++)
        {
            std::string result_name = name + "_result_" + std::to_string(i);
            all_results[i]->save(result_name);
        }
    }
    
    
    std::cout << "done decomposition of " << name << '\n';

    #if LOG_TIME == 1
    std::chrono::duration<float> duration_decomp = std::chrono::steady_clock::now() - start_decomp;
    std::cout << "[Duration] Decomposition: " << duration_decomp.count() << " s\n";
    time_o.open(Global_var::input_file + "_time.txt", std::ios::app);
    time_o << "Decomposition: " << duration_decomp.count() << " s\n";
    time_o.close();
    #endif  

    return 0;
}
