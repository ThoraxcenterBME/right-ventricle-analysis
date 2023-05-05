#include "strain_calculation.h"
#include "util.h"
#include <framework/mesh.h>
#include <framework/ray.h> 
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/reciprocal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <algorithm>
#include <functional>
#include <math.h>
#include <numeric>
#include <set>

// Calculates area strain
// Same formulation as: https://www.onlinejase.com/article/S0894-7317(14)00925-0/pdf
// Area strain defined as: https://doi.org/10.1093/ehjci/jeaa189
double area_strain(double ed_area, double es_area)
{
    return 100 * ((es_area - ed_area) / ed_area);
}

void set_regional_area_strain(Strain& strain)
{

    for (int i = 0; i < strain.es_areas.size(); i++) {
        strain.strain_values.push_back(area_strain(strain.ed_areas[i], strain.es_areas[i]));
    }
}

// Calculate the strain value of a vertex, along a specified axis 
double directional_strain_vertex(int vertex,
    std::vector<Vertex>& vertices_es,
    std::vector<Vertex>& vertices_ed,
    glm::vec3& axis)
{
    // Get current vertex from ED
    auto& current_vertex = vertices_ed[vertex];
    // Define origin as the current vertex's position
    // List of displacement vectors
    std::vector<double> u_vectors = {};
    auto length_axis = glm::length(axis) * glm::length(axis);

    for (int i = 0; i < current_vertex.ring.size(); i++) {
        // Find the displacement vectors of the opposite vertex in ES and ED.
        auto ed_vector = vertices_ed[current_vertex.ring[i]].position - vertices_ed[vertex].position;
        auto es_vector = vertices_es[current_vertex.ring[i]].position - vertices_es[vertex].position;

        // If too close to being perpendicular to the axis
        if (abs(glm::dot(ed_vector, axis) / length_axis)  < 1E-2) {
            continue;
        }

        // Map both onto the longitudinal axis and find length (nominal strain)
        auto proj_ed_vector = glm::length((glm::dot(ed_vector, axis) / length_axis) * axis);
        auto proj_es_vector = glm::length((glm::dot(es_vector, axis) / length_axis) * axis);

        // Strain percentage change
        auto strain = 100 * ((proj_es_vector - proj_ed_vector) / proj_ed_vector);

        // Subtract two vectors
        u_vectors.push_back(strain);
    }

    // Sum it up
    return std::accumulate(u_vectors.begin(), u_vectors.end(), 0.0) / static_cast<double>(current_vertex.ring.size());
}

// Returns regional longitudinal strain values
std::vector<double> set_regional_longitudinal_strain(Strain& strain, std::vector<Vertex>& vertices)
{
    std::map<Region, std::vector<double>> strain_regional = {};

    for (auto& v : vertices) {
        if (v.exclude)
            continue;

        strain_regional[v.region].push_back(v.long_strain);
    }

    std::vector<double> l_vals = {};
    for (auto const& r_long_strain : strain_regional) {
        // Calculate average longitudinal strain of that region
        l_vals.push_back(std::accumulate(r_long_strain.second.begin(), r_long_strain.second.end(), 0.0) / static_cast<double>(r_long_strain.second.size()));
    }
    return l_vals;
}

// Returns regional radial strain values
std::vector<double> set_regional_radial_strain(Strain& strain, std::vector<Vertex>& vertices)
{
    std::map<Region, std::vector<double>> strain_regional = {};

    for (auto& v : vertices) {
        if (v.exclude)
            continue;

        strain_regional[v.region].push_back(v.rad_strain);
    }

    std::vector<double> r_vals = {};
    for (auto const& r_long_strain : strain_regional) {
        // Calculate average longitudinal strain of that region
        r_vals.push_back(std::accumulate(r_long_strain.second.begin(), r_long_strain.second.end(), 0.0) / static_cast<double>(r_long_strain.second.size()));
    }
    return r_vals;
}

// Returns regional circumferential strain values
std::vector<double> set_regional_circumferential_strain(Strain& strain, std::vector<Vertex>& vertices)
{
    std::map<Region, std::vector<double>> strain_regional = {};

    for (auto& v : vertices) {
        if (v.exclude)
            continue;

        strain_regional[v.region].push_back(v.circ_strain);
    }

    std::vector<double> l_vals = {};
    for (auto const& r_long_strain : strain_regional) {
        // Calculate average longitudinal strain of that region
        l_vals.push_back(std::accumulate(r_long_strain.second.begin(), r_long_strain.second.end(), 0.0) / static_cast<double>(r_long_strain.second.size()));
    }
    return l_vals;
}

//bool include_in_strain(Vertex& vertex,
//    glm::vec3& l_axis,
//    std::vector<Ray>& normals,
//    std::map<int, std::vector<int>>& vertexToTri)
//{
//    // Check if any of the vertices in vertex ring have a normal perpendicular to long axis
//    for (auto const& i : vertex.ring) {
//        glm::vec3 n = vertex_normal(i, vertexToTri, normals);
//        float cosine = glm::dot(n, l_axis) / (glm::length(n) * glm::length(l_axis));
//        float ang = glm::acos(cosine);
//
//        // Do not include in calculations if perpendicular
//        if (std::fabs(ang - (M_PI / 2)) < 1E-2) {
//            return false;
//        }
//    }
//    return true;
//}

// Find longitudinal strain of the entire mesh
// Note: don't remove the `vertexToTri` argument 
double longitudinal_strain(std::vector<Vertex>& vertices_es,
    std::vector<Vertex>& vertices_ed,
    std::vector<glm::uvec3> triangles,
    std::map<int, std::vector<int>>& vertexToTri,
    glm::vec3& l_axis,
    Strain& strain)
{
    auto longitudinal_strain = 0.0;
    std::vector<Ray> normals = find_normals(vertices_ed, triangles);

    for (int i = 0; i < vertices_ed.size(); i++) {
        float strain_i = directional_strain_vertex(i, vertices_es, vertices_ed, l_axis);
        vertices_ed[i].long_strain = strain_i;
        longitudinal_strain += strain_i;
    }
    // Set the regional strain values
    auto reg_strain = set_regional_longitudinal_strain(strain, vertices_ed);
    strain.l_strain_values = reg_strain;

    // Return global strain values
    return longitudinal_strain / vertices_ed.size();
}

// Radial strain
double radial_strain(std::vector<Vertex>& vertices_es,
    std::vector<Vertex>& vertices_ed,
    glm::vec3& r_axis,
    Strain& strain)
{
    auto radial_strain = 0.0;

    for (int i = 0; i < vertices_ed.size(); i++) {
        float strain_i = directional_strain_vertex(i, vertices_es, vertices_ed, r_axis);
        vertices_ed[i].long_strain = strain_i;
        radial_strain += strain_i;
    }
    // Set the regional strain values
    auto reg_strain = set_regional_longitudinal_strain(strain, vertices_ed);
    strain.r_strain_values = reg_strain;

    // Return global strain values
    return radial_strain / vertices_ed.size();
}

double circumferential_strain(std::vector<Vertex>& vertices_es,
    std::vector<Vertex>& vertices_ed,
    glm::vec3& c_axis,
    Strain& strain)
{
    auto circumferential_strain = 0.0;

    for (int i = 0; i < vertices_ed.size(); i++) {
        float strain_i = directional_strain_vertex(i, vertices_es, vertices_ed, c_axis);
        vertices_ed[i].long_strain = strain_i;
        circumferential_strain += strain_i;
    }
    // Set the regional strain values
    auto reg_strain = set_regional_longitudinal_strain(strain, vertices_ed);
    strain.c_strain_values = reg_strain;

    // Return global strain values
    return circumferential_strain / vertices_ed.size();
}

// Determines the longitudinal axis to be used
glm::vec3 find_long_axis(std::vector<Vertex>& vertices, int lower_point, int c1, int c2)
{
    glm::vec3& apex = vertices[lower_point].position;

    // Take the center points of the inflow/outflow tracts
    glm::vec3& a = vertices[c1].position;
    glm::vec3& b = vertices[c2].position;

    // Define top point to be the midpoint
    glm::vec3 center = glm::vec3(0.5f * (a.x + b.x), 0.5f * (a.y + b.y), 0.5f * (a.z + b.z));

    return center - apex;
}


// Determines the radial axis to be used 
glm::vec3 find_radial_axis(std::vector<Vertex>& vertices, glm::vec3& l_axis) 
{
    std::vector<glm::vec3> positions;
    std::transform(std::begin(vertices), std::end(vertices),
        std::back_inserter(positions),
        [](const Vertex& v) { return v.position; });

    // Find the center of mass
    const glm::vec3 center = std::accumulate(std::begin(positions), std::end(positions), glm::vec3(0.0f)) / static_cast<float>(positions.size());

    // Use Gram-Schmidt process to find a vector orthogonal to longitudinal axis 
    glm::vec3 u = glm::normalize(l_axis); 
    glm::vec3 v = glm::vec3(1, 0, 0); 
    glm::vec3 p = glm::dot(v, u) * u; 
    glm::vec3 w = v - p; 
    glm::vec3 radial_dir = glm::normalize(w); 

    return center + radial_dir; 
}

// Determines the circumferential axis to be used 
glm::vec3 find_circumferential_axis(std::vector<Vertex>& vertices, glm::vec3& l_axis, glm::vec3& r_axis) 
{
    std::vector<glm::vec3> positions;
    std::transform(std::begin(vertices), std::end(vertices),
        std::back_inserter(positions),
        [](const Vertex& v) { return v.position; });

    // Find the center of mass
    const glm::vec3 center = std::accumulate(std::begin(positions), std::end(positions), glm::vec3(0.0f)) / static_cast<float>(positions.size());
    glm::vec3 circ_dir = glm::normalize(glm::cross(l_axis, r_axis)); 

    return center + circ_dir; 
}
