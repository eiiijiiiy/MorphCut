#pragma once

#include "util.hpp"
#include "global_variables.hpp"
#include "Component.hpp"
#include "Node.hpp"
#include "bnb.hpp"
#include "odas/nlopt.hpp"


std::pair<float, float> detect_symmetry(
    const std::vector<Point> &pts, 
    const std::vector<std::vector<size_t>> &triangles)
{
    Solvernlopt solver;
    std::cout << "loading mesh...\n";
    solver.load_prob(pts, triangles, 4, 1.0);
    std::cout << "loaded mesh!\n";
    solver.solve(0, 1000);
    solver.validate_best_sol();  

    float time_cost = std::round(solver.timecost*1000)/1000;
    float rmse = std::round(solver.valid_rmse*1000)/1000;
    std::cout << "Time cost: " << time_cost << "s, RMSE: " << rmse << '\n';
    float rho = solver.best_solution[0];
    float theta = solver.best_solution[1];

//     // visualization for test
#if VIS == 1
    Surface_mesh vis_mesh;
    std::vector<Vertex_index> vtxs;
    for (auto pt : pts)
        vtxs.push_back(vis_mesh.add_vertex(pt));
    for (auto f : triangles)
    {
        std::vector<Vertex_index> f_vtxs;
        for (auto vid : f)
            f_vtxs.push_back(vtxs[vid]);
        vis_mesh.add_face(f_vtxs);
    }

    Point p1(rho * std::cos(theta) - 100 * std::sin(theta), 
            rho * std::sin(theta) + 100 * std::cos(theta),
            -100);
    Point p2(rho * std::cos(theta) + 100 * std::sin(theta), 
            rho * std::sin(theta) - 100 * std::cos(theta),
            -100);
    Point p3(rho * std::cos(theta) - 100 * std::sin(theta), 
            rho * std::sin(theta) + 100 * std::cos(theta),
            100);
    Point p4(rho * std::cos(theta) + 100 * std::sin(theta), 
            rho * std::sin(theta) - 100 * std::cos(theta),
            100);
    Vertex_index v1 = vis_mesh.add_vertex(p1);
    Vertex_index v2 = vis_mesh.add_vertex(p2);
    Vertex_index v3 = vis_mesh.add_vertex(p3);
    Vertex_index v4 = vis_mesh.add_vertex(p4);

    vis_mesh.add_face(v1, v2, v3);
    vis_mesh.add_face(v2, v4, v3);
    CGAL::draw(vis_mesh, "vertical symmetry plane");
    // ------ visualization for test
#endif

    return std::make_pair(rho, theta);
}

