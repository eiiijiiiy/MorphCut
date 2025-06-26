#pragma once

#include "global_variables.hpp"
#include "util.hpp"


void raycasting(
    const std::vector<Point> &query_points,
    SDF_Result &results,
    const size_t num_ray_per_pt = 5)
{
    assert (results.rows() == query_points.size());

    #pragma omp parallel for
    for (size_t i = 0; i < query_points.size(); ++i) {
        const Point& point = query_points[i];
        Eigen::ArrayXi is_in = Eigen::ArrayXi::Zero(num_ray_per_pt);

        for (size_t j = 0; j < num_ray_per_pt; ++j) {
            // Generate a random direction for the ray
            float vertical_angle = std::rand()/double(RAND_MAX) * (M_PI /3) + M_PI /6;
            float horizontal_angle = std::rand()/double(RAND_MAX) * 2 * M_PI;
            Vector direction (
                cos(vertical_angle) * cos(horizontal_angle),
                cos(vertical_angle) * sin(horizontal_angle),
                sin(vertical_angle));
            Ray ray_query(point, point + direction);
            try{
                int num_inter = Global_var::tree.number_of_intersected_primitives(ray_query);
                is_in(j) = num_inter % 2 == 0 ? 0 : 1;
            }
            catch(CGAL::Precondition_exception e)
            {
                // TODO: delete these points?
                std::cerr << "precondition exception: " << e.what() << '\n';
            }
            results(i) = float(is_in.sum()) / num_ray_per_pt > 0.5 ? 1 : 0;
        }
    }
}

void sample_points_in_convex_cell3(
    const LCC_3 &lcc,
    Dart_const_descriptor &cell_dart,
    const size_t num_points,
    SDF_Points &query_points)
{
    //1. create tetrehedron (r1, r2, r3, r4) from the vertices of the cell (one point and the other three vertices that are not on any plane that same as a given point)
    //2. sample points in the tetrahedron: r = T*(lambda_1, lambda_2, lambda_3) + r_4
    // T = 
    // [x1-x4, x2-x4, x3-x4; y1-y4, y2-y4, y3-y4; z1-z4, z2-z4, z3-z4]
    std::vector<std::pair<std::vector<Eigen::Vector3f>, float>> tetrahedra;
    std::vector<Point> vtxs;
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::const_iterator
        v = lcc.one_dart_per_incident_cell<0, 3>(cell_dart).begin();
        v != lcc.one_dart_per_incident_cell<0, 3>(cell_dart).end(); ++v)
        vtxs.push_back(lcc.point(v));
    
    Triangulation T(vtxs.begin(), vtxs.end());
    float total_volume = 0;
    for (Finite_cells_iterator cit = T.finite_cells_begin();
        cit != T.finite_cells_end(); ++cit)
    {
        std::vector<Eigen::Vector3f> tetra;
        for (int i = 0; i < 4; i++)
        {
            Point pt = cit->vertex(i)->point();
            tetra.push_back(Eigen::Vector3f(
                static_cast<float>(CGAL::to_double(pt.x())),
                static_cast<float>(CGAL::to_double(pt.y())),
                static_cast<float>(CGAL::to_double(pt.z()))));
        }
        float volume = volume_of_tetrahedron(tetra);
        tetrahedra.push_back({tetra, volume});
        total_volume += volume;
    }
    std::cout << "#tetrahedra: " << tetrahedra.size() << '\n';

    float unit_volume = total_volume / num_points;
    // num_points = num_tetrahedra * num_points_per_tetrahedra
    size_t num_sampled_pts = 0;
    
    for (size_t i = 0; i < tetrahedra.size(); i++)
    {
        const auto [vtxs, volume] = tetrahedra[i];
        int num_pts_this_tetra = int(std::ceil(volume / unit_volume));
        int num_pins = int(std::ceil(std::pow(num_pts_this_tetra, 1.0/3.0)));
        // 1^3; 2^3, 3^3, 4^3, 5^3, 6^3, 7^3, 8^3, 9^3, 10^3 = 1000
        Eigen::MatrixXf weights = make_para_comb(4, num_pins);
        assert(weights.rows() >= num_pts_this_tetra);
        for (size_t j = 0; j < num_pts_this_tetra; j++)
        {
            Eigen::RowVector3f pt = Eigen::RowVector3f::Zero();
            for (size_t k = 0; k < 4; k++)
                pt = pt.array() + vtxs[k].transpose().array() * weights(j, k);
            query_points.row(num_sampled_pts) = pt;
            num_sampled_pts++;
            if (num_sampled_pts >= num_points)
                break;
        }
        if (num_sampled_pts >= num_points)
            break;
    }
}

float query_in_ratio_of_convex_cell3(
    std::string save_title,
    Dart_const_descriptor cell_dart)
{
    float volume = volume_of_cell3(Global_var::lcc, cell_dart);
    size_t num_points_of_volume = std::ceil(volume/Global_var::sampled_unit_volume);
    size_t num_points = std::max(Global_var::min_points, num_points_of_volume);
    num_points = std::min(num_points, Global_var::max_points);

    // std::vector<Eigen::Vector3f> vtxs;
    // for (LCC_3::One_dart_per_incident_cell_range<0, 3>::const_iterator
    //     v = Global_var::lcc.one_dart_per_incident_cell<0, 3>(cell_dart).begin();
    //     v != Global_var::lcc.one_dart_per_incident_cell<0, 3>(cell_dart).end(); ++v)
    //     vtxs.push_back(Eigen::Vector3f(
    //         static_cast<float>(CGAL::to_double(Global_var::lcc.point(v).x())),
    //         static_cast<float>(CGAL::to_double(Global_var::lcc.point(v).y())),
    //         static_cast<float>(CGAL::to_double(Global_var::lcc.point(v).z()))));

    // size_t num_vtx = vtxs.size();
    // std::cout << "num_points: " << num_points << ", num_vtx=" << num_vtx << '\n';
    // Eigen::Vector3f centroid = get_centroid(vtxs);
    // std::vector<float> dist_to_centroid(num_vtx);
    // float max_dist = 0;
    // for (size_t i = 0; i < num_vtx; i++)
    // {
    //     dist_to_centroid[i] = (vtxs[i] - centroid).norm();
    //     if (dist_to_centroid[i] > max_dist)
    //         max_dist = dist_to_centroid[i];
    // }
    // std::vector<float> vtx_weights(num_vtx);
    // for (size_t i = 0; i < num_vtx; i++)
    //     vtx_weights[i] = dist_to_centroid[i] / max_dist;

    // // Eigen::MatrixXf::NullaryExpr random(num_points, num_vtx, [&](){return dis(gen);});
    // // Eigen::MatrixXf weights = Eigen::MatrixXf::Zero(num_points, num_vtx).unaryExpr([&](float dummy){return dis(gen);});
    // Eigen::MatrixXf weights = Eigen::MatrixXf::Zero(num_points + num_vtx, num_vtx);
    // // Eigen::MatrixXf weights = random;
    // // weights = weights.array() + 1.0f;
    // // weights = weights.array() / (weights.rowwise().sum().replicate(1, num_vtx)).array();
    // for (size_t i = 0; i < num_vtx; i++)
    // {
    //     for (size_t j = 0; j < num_points/2; j++)
    //     {
    //         std::random_device rd;
    //         std::mt19937 gen(rd());  //here you could also set a seed
    //         float hi = (1 - weights.row(j).sum()) * vtx_weights[(j % num_vtx + i) % num_vtx];
    //         std::uniform_real_distribution<float> dis(0, hi);
    //         weights(j, (j % num_vtx + i) % num_vtx) = dis(gen);
    //     }
    // }

    // for (size_t i = 0; i < num_vtx; i++)
    // {
    //     for (size_t j = num_points/2; j < num_points; j++)
    //     {
    //         std::random_device rd;
    //         std::mt19937 gen(rd());  //here you could also set a seed
    //         std::uniform_real_distribution<float> dis(0, 1);
    //         weights(j, i) = dis(gen);
    //     }
    // }

    // SDF_Points query_points(num_points + num_vtx, 3);
    // for (size_t i = 0; i < num_points; i++)
    // {
    //     weights.row(i) = weights.row(i).array() / weights.row(i).sum();
    //     Eigen::RowVector3f pt = Eigen::RowVector3f::Zero();
    //     for (size_t j = 0; j < num_vtx; j++)
    //         pt = pt.array() + vtxs[j].transpose().array() * weights(i, j);
    //     query_points.row(i) = pt;
    // }
    // float other_weight = 0.1 / (num_vtx-1);
    // for (size_t i = 0; i < num_vtx; i++)
    // {
    //     Eigen::RowVector3f pt = Eigen::RowVector3f::Zero();
    //     for (size_t j = 0; j < num_vtx; j++)
    //         if (j != i)
    //             pt = pt.array() + vtxs[j].transpose().array() * other_weight;
    //         else
    //             pt = pt.array() + vtxs[j].transpose().array() * 0.9;
    //     query_points.row(num_points + i) = pt;
    // }

    SDF_Points query_points(num_points, 3);
    sample_points_in_convex_cell3(Global_var::lcc, cell_dart, num_points, query_points);

    // save points for test
    // std::string pt_path = Global_var::input_file + ".cut_" + save_title + ".xyz";
    // std::ofstream pt_o(pt_path);
    // for (size_t i = 0; i < num_points; i++)
    //     pt_o << query_points(i, 0) << ' ' << query_points(i, 1) << ' ' << query_points(i, 2) << '\n';
    // pt_o.close();
    SDF_Result results (query_points.rows());
    if (Global_var::use_raycasting)
    {
        std::vector<Point> query_points_cgal(query_points.rows());
        for (size_t i = 0; i < query_points.rows(); i++)
            query_points_cgal[i] = Point(
                query_points.row(i)[0], query_points.row(i)[1], query_points.row(i)[2]);
        
        raycasting(query_points_cgal, results);
    }
    else
    {
        Global_var::sdf.contains(query_points, results);
    }

    Global_var::cell_sdf[cell_dart] = std::make_pair(query_points, results);
    std::cout << "num in: " << results.array().count() << '\n';
    int num_in = results.array().count();
    float ratio = float(num_in) / (query_points.rows());
    return ratio;
}

Dart_descriptor box_cell()
{
    FT dx = Global_var::box.xmax() - Global_var::box.xmin();
    FT dy = Global_var::box.ymax() - Global_var::box.ymin();
    FT dz = Global_var::box.zmax() - Global_var::box.zmin();
    // FT radius = FT(0.5) * CGAL::sqrt(dx * dx + dy * dy + dz * dz);
    FT radius = FT(0.5) * (dx + dy + dz);
    // FT offset = radius * FT(0.01);
    FT offset = 0.1;

    // make the box larger to ensure all points are enclosed.
    FT xmin = Global_var::box.xmin() - offset, xmax = Global_var::box.xmax() + offset;
    FT ymin = Global_var::box.ymin() - offset, ymax = Global_var::box.ymax() + offset;
    FT zmin = Global_var::box.zmin() - offset, zmax = Global_var::box.zmax() + 2 * Global_var::cell_size;

    Dart_descriptor d = Global_var::lcc.make_hexahedron(
            LCC_3::Point(xmax, ymin, zmin),
            LCC_3::Point(xmax, ymax, zmin),
            LCC_3::Point(xmin, ymax, zmin),
            LCC_3::Point(xmin, ymin, zmin),
            LCC_3::Point(xmin, ymin, zmax),
            LCC_3::Point(xmax, ymin, zmax),
            LCC_3::Point(xmax, ymax, zmax),
            LCC_3::Point(xmin, ymax, zmax));

    Global_var::lcc.set_attribute<3>(d, Global_var::lcc.create_attribute<3>()); // set as default = 0
    Global_var::lcc.info<3>(d) = std::make_pair(0, false);

    // set the attribute of edges
    LCC_3::One_dart_per_incident_cell_range<1, 3> edges = Global_var::lcc.one_dart_per_incident_cell<1, 3>(d);
    int eid = 0;
    for (LCC_3::One_dart_per_incident_cell_range<1, 3>::iterator e = edges.begin(); e != edges.end(); ++e)
    {
        Global_var::lcc.set_attribute<1>(e, Global_var::lcc.create_attribute<1>());
        Global_var::lcc.info<1>(e) = eid;
        eid++;
    }

    // set the attribute of faces
    LCC_3::One_dart_per_incident_cell_range<2, 3> lcc_faces = Global_var::lcc.one_dart_per_incident_cell<2, 3>(d);
    int fid = 0;
    for (LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator f = lcc_faces.begin(); f != lcc_faces.end(); ++f)
    {
        Global_var::lcc.set_attribute<2>(f, Global_var::lcc.create_attribute<2>());
                
        const Plane* plane = normalize_plane(new Plane(
            Global_var::lcc.point(f), 
            Global_var::lcc.point(Global_var::lcc.previous(f)), 
            Global_var::lcc.point(Global_var::lcc.next(f))));

        Global_var::face_planes.insert({fid, plane});
        Global_var::lcc.info<2>(f) = std::make_pair(fid, plane);
        Global_var::all_face_planes_copy.insert(plane);
        fid++;
    }

    // set the attribute of vertices
    int vid = 0;
    LCC_3::One_dart_per_incident_cell_range<0, 3> vertices = Global_var::lcc.one_dart_per_incident_cell<0, 3>(d);
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator v = vertices.begin(); v != vertices.end(); ++v)
    {
        Global_var::lcc.info<0>(v) = vid;
        vid++;
    }

    return d;
}




// sample grid points
Dart_descriptor sample_grid_points(
    const float &cell_size)
{
    FT dx = Global_var::box.xmax() - Global_var::box.xmin();
    FT dy = Global_var::box.ymax() - Global_var::box.ymin();
    FT dz = Global_var::box.zmax() - Global_var::box.zmin();
    // FT radius = FT(0.5) * CGAL::sqrt(dx * dx + dy * dy + dz * dz);
    FT radius = FT(0.5) * (dx + dy + dz);
    // FT offset = radius * FT(0.01);
    FT offset = 0.1;

    // make the box larger to ensure all points are enclosed.
    FT xmin = Global_var::box.xmin() - offset, xmax = Global_var::box.xmax() + offset;
    FT ymin = Global_var::box.ymin() - offset, ymax = Global_var::box.ymax() + offset;
    FT zmin = Global_var::box.zmin() - offset, zmax = Global_var::box.zmax() + 2 * cell_size;

    Dart_descriptor d = Global_var::lcc.make_hexahedron(
            LCC_3::Point(xmax, ymin, zmin),
            LCC_3::Point(xmax, ymax, zmin),
            LCC_3::Point(xmin, ymax, zmin),
            LCC_3::Point(xmin, ymin, zmin),
            LCC_3::Point(xmin, ymin, zmax),
            LCC_3::Point(xmax, ymin, zmax),
            LCC_3::Point(xmax, ymax, zmax),
            LCC_3::Point(xmin, ymax, zmax));

    Global_var::lcc.set_attribute<3>(d, Global_var::lcc.create_attribute<3>()); // set as default = 0
    Global_var::lcc.info<3>(d) = std::make_pair(0, false);

    // set the attribute of edges
    LCC_3::One_dart_per_incident_cell_range<1, 3> edges = Global_var::lcc.one_dart_per_incident_cell<1, 3>(d);
    int eid = 0;
    for (LCC_3::One_dart_per_incident_cell_range<1, 3>::iterator e = edges.begin(); e != edges.end(); ++e)
    {
        Global_var::lcc.set_attribute<1>(e, Global_var::lcc.create_attribute<1>());
        Global_var::lcc.info<1>(e) = eid;
        eid++;
    }

    // set the attribute of faces
    LCC_3::One_dart_per_incident_cell_range<2, 3> lcc_faces = Global_var::lcc.one_dart_per_incident_cell<2, 3>(d);
    int fid = 0;
    for (LCC_3::One_dart_per_incident_cell_range<2, 3>::iterator f = lcc_faces.begin(); f != lcc_faces.end(); ++f)
    {
        Global_var::lcc.set_attribute<2>(f, Global_var::lcc.create_attribute<2>());
                
        const Plane* plane = normalize_plane(new Plane(
            Global_var::lcc.point(f), 
            Global_var::lcc.point(Global_var::lcc.previous(f)), 
            Global_var::lcc.point(Global_var::lcc.next(f))));

        Global_var::face_planes.insert({fid, plane});
        Global_var::lcc.info<2>(f) = std::make_pair(fid, plane);
        Global_var::all_face_planes_copy.insert(plane);
        fid++;
    }

    // set the attribute of vertices
    int vid = 0;
    LCC_3::One_dart_per_incident_cell_range<0, 3> vertices = Global_var::lcc.one_dart_per_incident_cell<0, 3>(d);
    for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator v = vertices.begin(); v != vertices.end(); ++v)
    {
        Global_var::lcc.info<0>(v) = vid;
        vid++;
    }

    int num_x = static_cast<int>(CGAL::to_double((xmax - xmin) / cell_size)) + 1;
    int num_y = static_cast<int>(CGAL::to_double((ymax - ymin) / cell_size)) + 1;
    int num_z = static_cast<int>(CGAL::to_double((zmax - zmin) / cell_size)) + 1;

    Global_var::grid_points.clear();
    // Global_var::grid_points.resize(total_sampled_grid_pts);

    for (size_t xi = 0; xi < num_x; ++xi) {
        for (size_t yi = 0; yi < num_y; ++yi) {
            for (size_t zi = 0; zi < num_z; ++zi) {
                Global_var::grid_points.push_back(Point(
                    xmin + FT(xi) * cell_size,
                    ymin + FT(yi) * cell_size,
                    zmin + FT(zi) * cell_size));
            }
        }
    }
    int sparse_grid_pts = Global_var::grid_points.size();
    std::cout << "sampled " << sparse_grid_pts << " grid points\n";

    // denser sampling around the planes
    float smaller_cell_size = cell_size * 0.25;
    std::cout << "smaller cell size: " << smaller_cell_size << '\n';
    num_x = static_cast<int>(CGAL::to_double((xmax - xmin) / smaller_cell_size)) + 1;
    num_y = static_cast<int>(CGAL::to_double((ymax - ymin) / smaller_cell_size)) + 1;
    num_z = static_cast<int>(CGAL::to_double((zmax - zmin) / smaller_cell_size)) + 1;

    std::vector<bool> near_plane(num_x * num_y * num_z, false);
    // #pragma omp parallel for num_threads(16)
    // size_t count_plane = 0;
    // for (auto &[plane, region] : Global_var::plane_regions)
    // {
    //     count_plane++;
    //     std::cout << "\rplane " << count_plane << '/' << Global_var::plane_regions.size() << std::flush;
        
    //     for (auto f : region)
    //     {
    //         std::vector<Point> pts;
    //         Halfedge_index e = Global_var::input_mesh.halfedge(f);
    //         Halfedge_index begin_e = e;
    //         do {
    //             Vertex_index v = Global_var::input_mesh.source(e);
    //             pts.push_back(Global_var::input_mesh.point(v));
    //             e = Global_var::input_mesh.next(e);
    //         } while (e != begin_e);
    //         auto bbox = CGAL::bounding_box(pts.begin(), pts.end());
    //         FT xmin_f = bbox.xmin() - smaller_cell_size, xmax_f = bbox.xmax() + smaller_cell_size;
    //         FT ymin_f = bbox.ymin() - smaller_cell_size, ymax_f = bbox.ymax() + smaller_cell_size;
    //         FT zmin_f = bbox.zmin() - smaller_cell_size, zmax_f = bbox.zmax() + smaller_cell_size;

    //         int num_x_f = static_cast<int>(CGAL::to_double((xmax_f - xmin_f) / smaller_cell_size)) + 1;
    //         int num_y_f = static_cast<int>(CGAL::to_double((ymax_f - ymin_f) / smaller_cell_size)) + 1;
    //         int num_z_f = static_cast<int>(CGAL::to_double((zmax_f - zmin_f) / smaller_cell_size)) + 1;
    //         std::cout << "num_x_f: " << num_x_f << " num_y_f: " << num_y_f << " num_z_f: " << num_z_f << '\n';
    //         for (size_t xi = 0; xi < num_x_f; ++xi) {
    //             for (size_t yi = 0; yi < num_y_f; ++yi) {
    //                 for (size_t zi = 0; zi < num_z_f; ++zi) {
    //                     Point pt (xmin_f + FT(xi) * smaller_cell_size,
    //                         ymin_f + FT(yi) * smaller_cell_size,
    //                         zmin_f + FT(zi) * smaller_cell_size);
    //                     Global_var::grid_points.push_back(pt);
    //                 }
    //             }
    //         }
    //     }
    // }
    // Surface_mesh mesh;
    // if (Global_var::use_repair) 
    //     mesh = Global_var::repaired_mesh;
    // else 
    //     mesh = Global_var::input_mesh;

    // for (Face_index f : faces(mesh))
    // {
    //     std::vector<Point> pts;
    //     Halfedge_index e = mesh.halfedge(f);
    //     Halfedge_index begin_e = e;
    //     do {
    //         Vertex_index v = mesh.source(e);
    //         pts.push_back(mesh.point(v));
    //         e = mesh.next(e);
    //     } while (e != begin_e);
    //     auto bbox = CGAL::bounding_box(pts.begin(), pts.end());
    //     FT xmin_f = bbox.xmin() - smaller_cell_size, xmax_f = bbox.xmax() + smaller_cell_size;
    //     FT ymin_f = bbox.ymin() - smaller_cell_size, ymax_f = bbox.ymax() + smaller_cell_size;
    //     FT zmin_f = bbox.zmin() - smaller_cell_size, zmax_f = bbox.zmax() + smaller_cell_size;

    //     int num_x_f = static_cast<int>(CGAL::to_double((xmax_f - xmin_f) / smaller_cell_size)) + 1;
    //     int num_y_f = static_cast<int>(CGAL::to_double((ymax_f - ymin_f) / smaller_cell_size)) + 1;
    //     int num_z_f = static_cast<int>(CGAL::to_double((zmax_f - zmin_f) / smaller_cell_size)) + 1;
    //     std::cout << "num_x_f: " << num_x_f << " num_y_f: " << num_y_f << " num_z_f: " << num_z_f << '\n';
    //     for (size_t xi = 0; xi < num_x_f; ++xi) {
    //         for (size_t yi = 0; yi < num_y_f; ++yi) {
    //             for (size_t zi = 0; zi < num_z_f; ++zi) {
    //                 Point pt (xmin_f + FT(xi) * smaller_cell_size,
    //                     ymin_f + FT(yi) * smaller_cell_size,
    //                     zmin_f + FT(zi) * smaller_cell_size);
    //                 Global_var::grid_points.push_back(pt);
    //             }
    //         }
    //     }
    // }
    // #pragma omp parallel for
    // for (size_t xi = 0; xi < num_x; ++xi) {
    //     for (size_t yi = 0; yi < num_y; ++yi) {
    //         for (size_t zi = 0; zi < num_z; ++zi) {
    //             Point pt (xmin + FT(xi) * smaller_cell_size,
    //                 ymin + FT(yi) * smaller_cell_size,
    //                 zmin + FT(zi) * smaller_cell_size);
    //             for (auto plane : Global_var::rg_planes)
    //             {
    //                 if (dist_pt_to_plane(pt, *plane) < cell_size*0.5)
    //                 {
    //                     near_plane[xi * num_y * num_z + yi * num_z + zi] = true;
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    // }

    // for (size_t xi = 0; xi < num_x; ++xi) {
    //     for (size_t yi = 0; yi < num_y; ++yi) {
    //         for (size_t zi = 0; zi < num_z; ++zi) {
    //             if (near_plane[xi * num_y * num_z + yi * num_z + zi])
    //             {
    //                 Global_var::grid_points.push_back(Point(
    //                     xmin + FT(xi) * smaller_cell_size,
    //                     ymin + FT(yi) * smaller_cell_size,
    //                     zmin + FT(zi) * smaller_cell_size));
    //             }
    //         }
    //     }
    // }

    // std::cout << "sampled " << Global_var::grid_points.size() - sparse_grid_pts 
    //         << "denser grid points around planes\n";
    return d;
}

void vote_in_out(
    const size_t num_ray_per_pt,
    const float in_thresh)
{
    // constructs AABB tree
    Tree tree;
    std::cout << "use repair: " << Global_var::use_repair << "\n";
    if (Global_var::use_repair)
        tree = Tree(faces(Global_var::repaired_mesh).first, 
                    faces(Global_var::repaired_mesh).second, Global_var::repaired_mesh);
    else
        tree = Tree(faces(Global_var::input_mesh).first, 
                    faces(Global_var::input_mesh).second, Global_var::input_mesh);

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
    sdf::SDF sdf(pts, tris, true);
    
    SDF_Points grid_pts(Global_var::grid_points.size(), 3);
    for (size_t i = 0; i < Global_var::grid_points.size(); ++i) {
        Eigen::RowVector3f pt {{
            static_cast<float>(CGAL::to_double(Global_var::grid_points[i].x())),
            static_cast<float>(CGAL::to_double(Global_var::grid_points[i].y())),
            static_cast<float>(CGAL::to_double(Global_var::grid_points[i].z()))}};
        grid_pts.row(i) = pt;
    }
    SDF_Result results (Global_var::grid_points.size());
    sdf.contains(grid_pts, results);
    Global_var::vote_map = results.array();
    // // // TODO: how to deal with the bottom face issue?

    // // ray intersections
    // // Eigen::initParallel();
    // std::vector<int> vote_map_(Global_var::grid_points.size(), 0);
    // std::cout << "begin checking the in/out of sampled points\n";

    // #pragma omp parallel for
    // for (size_t i = 0; i < Global_var::grid_points.size(); ++i) {
    //     // std::cout << "\rgrid point " << i << '/' << Global_var::grid_points.size() << std::flush;
    //     const Point& grid_point = Global_var::grid_points[i];
    //     Eigen::ArrayXi is_in = Eigen::ArrayXi::Zero(num_ray_per_pt);

    //     for (size_t j = 0; j < num_ray_per_pt; ++j) {
    //         // Generate a random direction for the ray
    //         float vertical_angle = std::rand()/double(RAND_MAX) * (M_PI /3) + M_PI /6;
    //         float horizontal_angle = std::rand()/double(RAND_MAX) * 2 * M_PI;
    //         Vector direction (
    //             cos(vertical_angle) * cos(horizontal_angle),
    //             cos(vertical_angle) * sin(horizontal_angle),
    //             sin(vertical_angle));
    //         Ray ray_query(grid_point, grid_point + direction);
    //         try{
    //             is_in(j) = tree.number_of_intersected_primitives(ray_query) % 2 == 0 ? 0 : 1;
    //         }
    //         catch(CGAL::Precondition_exception e)
    //         {
    //             // TODO: delete these points?
    //             std::cerr << "precondition exception: " << e.what() << '\n';
    //         }
            
    //     }
    //     vote_map_[i] = float(is_in.sum()) / num_ray_per_pt > in_thresh ? 1 : 0;
    //     // Global_var::vote_map(i) = float(is_in.sum()) / num_ray_per_pt > in_thresh ? 1 : 0;
    // }

    // // copy vote_map_ to Global_var::vote_map
    // Global_var::vote_map = Eigen::Map<Eigen::ArrayXi>(vote_map_.data(), vote_map_.size());
}

void pre_locate_grid_points_of_planes(
    const std::vector<const Plane*> &planes,
    std::unordered_map<const Plane*, Eigen::ArrayXi> &on_plane_pos)
{
    size_t num_grid_points = Global_var::grid_points.size();
    for (const auto &plane : planes) {
        Eigen::ArrayXi in_out = Eigen::ArrayXi::Zero(num_grid_points);
        for (size_t i = 0; i < num_grid_points; ++i) {
            const Point& grid_point = Global_var::grid_points[i];
            in_out(i) = plane->oriented_side(grid_point) == CGAL::ON_POSITIVE_SIDE ? 1 : 0;
        }
        on_plane_pos[plane] = in_out;
    }
}

void pre_locate_grid_points_of_planes()
{
    pre_locate_grid_points_of_planes(Global_var::rg_planes, Global_var::on_rg_plane_pos);
}

std::pair<int, int> get_num_in_out_of_cell(
    const std::vector<int> &in_cell_pts,
    const Eigen::Array<bool, Eigen::Dynamic, 1> &vote_map)
{
    // assert (in_cell_pts.size() > 0);
    int num_in = vote_map(in_cell_pts).count();
    return std::make_pair(num_in, in_cell_pts.size() - num_in);
}

float get_in_out_ratio_of_cell(
    const int num_in_pts,
    Dart_const_descriptor cell)
{
    float cell_volume = volume_of_cell3(Global_var::lcc, cell);
    float in_volume = num_in_pts * Global_var::sampled_unit_volume;
    float ratio = in_volume / cell_volume;
    ratio = std::clamp(ratio, 0.0f, 1.0f);
    return ratio;
}

std::pair<float, float> get_in_out_volume_of_cell(
    const int num_in_pts,
    Dart_const_descriptor cell)
{
    float cell_volume = volume_of_cell3(Global_var::lcc, cell);
    float in_volume = num_in_pts * Global_var::sampled_unit_volume;
    if (in_volume < cell_volume) return std::make_pair(in_volume, cell_volume - in_volume);
    else return std::make_pair(cell_volume, 0);
}

// float compute_plane_f1(
//     const Eva_cut& delta_cut)
// {
//     const auto &[delta_in_int, delta_int] = delta_cut;
//     const auto &[current_in_int, current_int] = Global_var::current_cut;
//     int updated_in_int = current_in_int + delta_in_int; // TP
//     int updated_int = current_int + delta_int; // P 
//     int update_out_int = updated_int - updated_in_int; // FP
//     int update_in_ext = Global_var::total_in - updated_in_int; // FN

//     float precision = float(updated_in_int) / (updated_in_int + update_out_int);
//     float recall = float(updated_in_int) / (updated_in_int + update_in_ext);
//     float F1 = 2 * precision * recall / (precision + recall);

//     return F1;
// }

// float compute_plane_delta_f1(
//     const Eva_cut& delta_cut)
// {
//     float f1_before = compute_plane_f1(std::make_pair(0, 0));
//     float f1_after = compute_plane_f1(delta_cut);
//     return f1_after - f1_before;
// }


// void print_precision_recall(
//     const Eva_cut& delta_cut)
// {
//     const auto &[delta_in_int, delta_int] = delta_cut;
//     const auto &[current_in_int, current_int] = Global_var::current_cut;
//     int updated_in_int = current_in_int + delta_in_int; // TP
//     int updated_int = current_int + delta_int; // P 
//     int update_out_int = updated_int - updated_in_int; // FP
//     int update_in_ext = Global_var::total_in - updated_in_int; // FN

//     float precision = float(updated_in_int) / (updated_in_int + update_out_int);
//     float recall = float(updated_in_int) / (updated_in_int + update_in_ext);
    
//     std::cout << "TP: " << updated_in_int << " FP: " << update_out_int << " FN: " << update_in_ext << '\n'
//               << " precision: " << precision << " recall: " << recall << '\n';
// }


// float compute_cell_score(
//     const int num_in_pts,
//     Dart_const_descriptor cell_dart,
//     const bool is_in)
// {
//     const auto [volume_in, volume_out] = get_in_out_volume_of_cell(num_in_pts, cell_dart);
//     if (is_in)
//         return volume_out;
//     else
//         return volume_in;
// }

float compute_cell_score(
    Dart_const_descriptor cell_dart,
    const float in_ratio,
    const bool is_in)
{
    float volume = volume_of_cell3(Global_var::lcc, cell_dart);
    if (is_in)
        return volume * (1 - in_ratio);
    else
        return volume * in_ratio;
}


void iterative_cut_for_reconstruction()
{
    Dart_const_descriptor const_cell_dart;
    // 1. find the top planes of all the cells
    std::vector<std::pair<Dart_descriptor, float>> cells_to_cut;
    std::vector<Dart_descriptor> remove_cells;
    float total_volume = 0;
    for (auto &[cell, to_cut_planes_of_the_cell] : Global_var::to_cut_planes)
    {
        total_volume += volume_of_cell3(Global_var::lcc, cell);
        if (to_cut_planes_of_the_cell.size() > 0)
        {
            std::sort(to_cut_planes_of_the_cell.begin(), to_cut_planes_of_the_cell.end(), 
                [](const std::pair<const Plane*, float> &a, const std::pair<const Plane*, float> &b) -> bool
                {
                    return a.second > b.second;
                });
            std::string pt_title = "before_cut"; 
            float ratio = query_in_ratio_of_convex_cell3(pt_title, cell);
            Global_var::lcc.info<3>(cell).second = ratio > 0.5;
            float score = compute_cell_score(cell, ratio, Global_var::lcc.info<3>(cell).second);
            cells_to_cut.push_back(std::make_pair(cell, score));
        }
        else
            remove_cells.push_back(cell);
    }
    std::cout << "total volume: " << total_volume << '\n';
    

    for (auto cell : remove_cells)
        Global_var::to_cut_planes.erase(cell);
        
    std::cout << "initial to cut cells: " << cells_to_cut.size() << '\n';

    // 2. sort the cells
    std::sort(cells_to_cut.begin(), cells_to_cut.end(),
        [](auto &a, auto &b) {
            return a.second > b.second;  // compare the score of cell
        });
    
    #if DEBUG_RECONSTRUCTION == 1
    std::cout << "# planes to cut: " << top_planes_of_cells.size() << '\n';
    for (const auto &[cell, eva] : top_planes_of_cells)
        std::cout << "cell score: " << cell_scores[cell] << " plane score: " << eva.first << '\n';
    #endif
    
    // draw_sorting(top_planes_of_cells, Global_var::to_cut_planes, Global_var::lcc, Global_var::input_mesh, Global_var::plane_regions);
    // sort the cells by the score of the top plane
    size_t num_cut = 0;
    std::string in_sdf_file, out_sdf_file;
    while (cells_to_cut.size() > 0)
    {
        // draw_sorting(top_planes_of_cells, Global_var::to_cut_planes, Global_var::lcc, Global_var::input_mesh, Global_var::plane_regions);
        auto [top_cell_to_cut, score] = cells_to_cut[0];
        std::cout << "* #to cut cells=" << cells_to_cut.size() << "(top score=" << score << ")\n";
        std::vector<std::pair<const Plane*, float>> &to_cut_planes_of_this_cell = Global_var::to_cut_planes[top_cell_to_cut];
        const Plane* cut_plane = to_cut_planes_of_this_cell[0].first;

        std::vector<Dart_descriptor> pos_cells, neg_cells;
        int result = cut_cell3_by_plane(Global_var::lcc, top_cell_to_cut, cut_plane, pos_cells, neg_cells, true);
        if (result == EXIT_FAILURE)
        {
            std::cout << "cut failed\n";
            std::cout << "to_cut_planes_of_this_cell.size(): " << to_cut_planes_of_this_cell.size() << '\n';
            to_cut_planes_of_this_cell.erase(to_cut_planes_of_this_cell.begin());
            if (to_cut_planes_of_this_cell.size() == 0)
            {
                std::cout << "no plane to cut\n";
                Global_var::to_cut_planes.erase(top_cell_to_cut);
                cells_to_cut.erase(cells_to_cut.begin());
            }
            continue;
        }
        assert (pos_cells.size() == 1 && neg_cells.size() == 1);
        Global_var::all_face_planes_copy.insert(cut_plane);
        num_cut++;
        std::cout << "cut " << num_cut << " times\n";
        cells_to_cut.erase(cells_to_cut.begin());
        Global_var::cell_sdf.erase(top_cell_to_cut);

        // const auto [pos_d, neg_d] = cut_cells;
        const auto pos_d = pos_cells[0], neg_d = neg_cells[0];
        std::string pos_title = "pos_cell_" + std::to_string(num_cut);
        std::string neg_title = "neg_cell_" + std::to_string(num_cut);
        float in_ratio_pos = query_in_ratio_of_convex_cell3(pos_title, pos_d),
              in_ratio_neg = query_in_ratio_of_convex_cell3(neg_title, neg_d);
        Global_var::lcc.info<3>(pos_d).second = in_ratio_pos > 0.5;
        Global_var::lcc.info<3>(neg_d).second = in_ratio_neg > 0.5;
        float score_pos = compute_cell_score(pos_d, in_ratio_pos, Global_var::lcc.info<3>(pos_d).second);
        float score_neg = compute_cell_score(neg_d, in_ratio_neg, Global_var::lcc.info<3>(neg_d).second);
        std::cout << "in_ratio_pos=" << in_ratio_pos << ", score_pos=" << score_pos << '\n'
                    << "in_ratio_neg=" << in_ratio_neg << ", score_neg=" << score_neg << '\n';

        // std::string cut_file_ = Global_var::input_file + ".cur_cut.off";
        // write_lcc_as_off(Global_var::lcc, cut_file_, true);
        // std::string cut_file_all_lcc_ = Global_var::input_file + ".cur_cut_all.off";
        // write_lcc_as_off(Global_var::lcc, cut_file_all_lcc_, false);

        // in_sdf_file = Global_var::input_file + ".cut_" + std::to_string(num_cut) + "_in_sdf.xyz";
        // out_sdf_file = Global_var::input_file + ".cut_" + std::to_string(num_cut) + "_out_sdf.xyz";
        // std::ofstream in_sdf_f(in_sdf_file), out_sdf_f(out_sdf_file);
        // for (const auto& [cell, query] : Global_var::cell_sdf)
        // {
        //     const auto &[query_pts, query_results] = query;
        //     for (size_t i = 0; i < query_pts.rows(); ++i)
        //     {
        //         if (query_results(i))
        //             in_sdf_f << query_pts.row(i) << '\n';
        //         else
        //             out_sdf_f << query_pts.row(i) << '\n';
        //     }
        // }
        // in_sdf_f.close();
        // out_sdf_f.close();
        
        // sort the top planes of the pos cells
        to_cut_planes_of_this_cell.erase(to_cut_planes_of_this_cell.begin());
        std::vector<std::pair<const Plane*, float>> to_cut_planes_for_pos, to_cut_planes_for_neg;
        // std::cout << "num_cut: " << num_cut << '\n';
        // std::cout << "to_cut_planes_of_this_cell.size(): " << to_cut_planes_of_this_cell.size() << '\n';
        for (const auto &[plane, area] : to_cut_planes_of_this_cell)
        {   
            // bool ch_pt_in_pos = false, ch_pt_in_neg = false;
            // for (auto ch_pt : Global_var::plane_extended_convex_hulls.at(plane))
            // {
            //     if (cut_plane->oriented_side(ch_pt) == CGAL::ON_POSITIVE_SIDE)
            //         ch_pt_in_pos = true;
            //     else
            //         ch_pt_in_neg = true;
            // }
            // if (ch_pt_in_pos)
            //     to_cut_planes_for_pos.push_back(std::make_pair(plane, 1.0));
            // if (ch_pt_in_neg)
            //     to_cut_planes_for_neg.push_back(std::make_pair(plane, 1.0));
            if(!is_cell3_on_one_side(Global_var::lcc, pos_d, plane, Global_var::input_mesh, Global_var::plane_regions[plane],
                Global_var::dist_thresh, 3 * Global_var::dist_thresh))
                to_cut_planes_for_pos.push_back(std::make_pair(plane, area));
            if(!is_cell3_on_one_side(Global_var::lcc, neg_d, plane, Global_var::input_mesh, Global_var::plane_regions[plane],
                Global_var::dist_thresh, 3 * Global_var::dist_thresh))
                to_cut_planes_for_neg.push_back(std::make_pair(plane, area));
        }
        std::cout << "# to cut plane of pos: " << to_cut_planes_for_pos.size() << ", of neg: " << to_cut_planes_for_neg.size() << '\n';
        
        float volume_pos = volume_of_cell3(Global_var::lcc, pos_d), 
              volume_neg = volume_of_cell3(Global_var::lcc, neg_d);
        if (to_cut_planes_for_pos.size() > 0 && 
            (score_pos > 1e-3 || (in_ratio_pos <= 0.5 && volume_pos / total_volume > 0.05)))
        {
            std::sort(to_cut_planes_for_pos.begin(), to_cut_planes_for_pos.end(), 
                [](const std::pair<const Plane*, float> &a, const std::pair<const Plane*, float> &b) -> bool
                {
                    return a.second > b.second;
                });
            Global_var::to_cut_planes[pos_d] = to_cut_planes_for_pos;
            cells_to_cut.push_back(std::make_pair(pos_d, score_pos));
        }
        else
        {
            std::cout << "* No plane to cut POS_CELL (in: " << Global_var::lcc.info<3>(pos_d).second << '\n';
            if (score_pos > 1e-3)
                std::cout << "[warning!] mixing pos cell (score=" << score_pos <<  ")without planes to further cut\n";
        }
            
        if (to_cut_planes_for_neg.size() > 0 && (score_neg > 1e-3 || (in_ratio_neg <= 0.5 && volume_neg / total_volume > 0.05)))
        // if (to_cut_planes_for_neg.size() > 0 && ratio_neg_in > 1e-3 && ratio_neg_in < 1-1e-3)
        {
            std::sort(to_cut_planes_for_neg.begin(), to_cut_planes_for_neg.end(), 
                [](const std::pair<const Plane*, float> &a, const std::pair<const Plane*, float> &b) -> bool
                {
                    return a.second > b.second;
                });

            Global_var::to_cut_planes[neg_d] = to_cut_planes_for_neg;
            cells_to_cut.push_back(std::make_pair(neg_d, score_neg));
        }
        else
        {
            std::cout << "* No plane to cut NEG_CELL (in: " << Global_var::lcc.info<3>(neg_d).second << '\n';
            if (score_neg > 1e-3)
                std::cout << "[warning!] mixing neg cell (score=" << score_neg <<  ")without planes to further cut\n";
        }
            
        // total num of xx change, the score and order also change
        std::sort(cells_to_cut.begin(), cells_to_cut.end(),
            [](auto &a, auto &b) {
                return a.second > b.second;  // compare the score of cell
            });
        #if DEBUG_RECONSTRUCTION == 1
        std::cout << "# planes to cut: " << top_planes_of_cells.size() << '\n';
        for (const auto &[cell, eva] : top_planes_of_cells)
            std::cout << "cell score: " << cell_scores[cell] << " plane score: " << eva.first << '\n';
        #endif
    }
    // refine
    // std::cout << "refining ... \n";
    // for (LCC_3::One_dart_per_cell_range<3>::iterator
    //     it=Global_var::lcc.one_dart_per_cell<3>().begin(),
    //     itend=Global_var::lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    // {
    //     float ratio = query_in_ratio_of_convex_cell3(it, 1000);
    //     Global_var::lcc.info<3>(it).second = ratio > 0.5;
    // }

    std::cout << "iterative cut finished\n";

    // std::string cut_file = Global_var::input_file + ".cut_" + std::to_string(num_cut) + ".off";
    // write_lcc_as_off(Global_var::lcc, cut_file, true);
    // std::string cut_file_all_lcc = Global_var::input_file + ".cut_" + std::to_string(num_cut) + "_all_lcc.off";
    // write_lcc_as_off(Global_var::lcc, cut_file_all_lcc, false);
    
    std::cout << "clean face inside ... \n";
    std::vector<Dart_descriptor> faces_to_remove;
    size_t num_volume_before_clean = Global_var::lcc.number_of_attributes<3>();

    for (LCC_3::One_dart_per_cell_range<2>::iterator
        it=Global_var::lcc.one_dart_per_cell<2>().begin(),
        itend=Global_var::lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
    {
        if (Global_var::lcc.is_free<3>(it)) continue;
        Dart_descriptor oppo_face = Global_var::lcc.beta<3>(it);
        if (oppo_face == it) continue;
        if (Global_var::lcc.info<3>(it).second == Global_var::lcc.info<3>(oppo_face).second)
            faces_to_remove.push_back(it);
    }

    // Draw_volume_functor df(amark);
    // CGAL::draw(Global_var::lcc, "iterative cut", false, df);
    for (auto face : faces_to_remove)
        Global_var::lcc.remove_cell<2>(face); // automatic update volume index by onmerge functor
        
    // update volume
    float max_volume = 0;
    for(LCC_3::One_dart_per_cell_range<3>::iterator
        it= Global_var::lcc.one_dart_per_cell<3>().begin(),
        itend= Global_var::lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
        if (!Global_var::lcc.info<3>(it).second)
            Global_var::lcc.remove_cell<3>(it);
        else
        {
            Dart_const_descriptor in_cell = it;
            float v_in_cell = volume_of_cell3(Global_var::lcc, in_cell);
            if (v_in_cell > max_volume)
                max_volume = v_in_cell;
            else
                Global_var::lcc.remove_cell<3>(it);
        }
    }

    int num_volume_after_clean = Global_var::lcc.number_of_attributes<3>();
    std::cout << "number of volume before clean: " << num_volume_before_clean 
              << " after clean: " << num_volume_after_clean << '\n';
    assert(num_volume_before_clean > num_volume_after_clean);
    
    // update vertex, edge, face, and volume index
    Vertex_idx vtx_id = 0;
    for(LCC_3::One_dart_per_cell_range<0>::iterator
        it= Global_var::lcc.one_dart_per_cell<0>().begin(),
        itend= Global_var::lcc.one_dart_per_cell<0>().end(); it!=itend; ++it)
        Global_var::lcc.info<0>(it) = vtx_id++;
    
    Edge_idx eid = 0;
    for(LCC_3::One_dart_per_cell_range<1>::iterator
        it= Global_var::lcc.one_dart_per_cell<1>().begin(),
        itend= Global_var::lcc.one_dart_per_cell<1>().end(); it!=itend; ++it)
        Global_var::lcc.info<1>(it) = eid++;

    float surface_area = 0;
    Face_idx fid = 0;
    for(LCC_3::One_dart_per_cell_range<2>::iterator
        it= Global_var::lcc.one_dart_per_cell<2>().begin(),
        itend= Global_var::lcc.one_dart_per_cell<2>().end(); it!=itend; ++it)
    {
        Global_var::lcc.info<2>(it).first = fid++;
        surface_area += area_of_cell2(Global_var::lcc, it);
    }
    Global_var::surface_area = surface_area;
    std::cout << "surface area: " << surface_area << '\n';

    Volume_idx vid = 0;
    for(LCC_3::One_dart_per_cell_range<3>::iterator
        it= Global_var::lcc.one_dart_per_cell<3>().begin(),
        itend= Global_var::lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
            Global_var::lcc.info<3>(it).first = vid++;
    assert (vid == 1);
    Dart_const_descriptor singleton_cell = Global_var::lcc.one_dart_per_cell<3>().begin();
    total_volume = volume_of_cell3(Global_var::lcc, singleton_cell);
    std::cout << "total volume: " << total_volume << '\n';
    #if VIS == 1
    CGAL::draw(Global_var::lcc, "clean reconstruction");
    #endif

    std::string clean_off = Global_var::input_file + "_reconstruction.off";
    std::string clean_sta = Global_var::input_file + "_reconstruction.txt";
    write_lcc_as_off(Global_var::lcc, clean_off, true, true, clean_sta);

    // save cell sdf
    // in_sdf_file = Global_var::input_file + "_in_sdf.xyz";
    // out_sdf_file = Global_var::input_file + "_out_sdf.xyz";
    // std::ofstream in_sdf_f(in_sdf_file), out_sdf_f(out_sdf_file);
    // for (const auto& [cell, query] : Global_var::cell_sdf)
    // {
    //     const auto &[query_pts, query_results] = query;
    //     for (size_t i = 0; i < query_pts.rows(); ++i)
    //     {
    //         if (query_results(i))
    //             in_sdf_f << query_pts.row(i) << '\n';
    //         else
    //             out_sdf_f << query_pts.row(i) << '\n';
    //     }
    // }
    // in_sdf_f.close();
    // out_sdf_f.close();
}
