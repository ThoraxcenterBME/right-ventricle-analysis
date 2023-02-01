#include "util.h"
#include <framework/mesh.h>
#include <framework/ray.h> 
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/reciprocal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <math.h>
#include <numeric>
constexpr auto EPS = 1e-7;
constexpr auto M_PI = 3.14159265358979323846; 
 
/*
 * Used the source: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 * to find the volume of a 3D mesh.
 */
double find_volume(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices)
{
    double volume = 0.0;

    for (int i = 0; i < triangles.size(); i++) {
        glm::uvec3 triangle = triangles[i];

        Vertex& p1 = vertices[triangle[0]];
        Vertex& p2 = vertices[triangle[1]];
        Vertex& p3 = vertices[triangle[2]];

        volume += (1.0 / 6.0) * (-p3.position[0] * p2.position[1] * p1.position[2] + p2.position[0] * p3.position[1] * p1.position[2] + p3.position[0] * p1.position[1] * p2.position[2] - p1.position[0] * p3.position[1] * p2.position[2] - p2.position[0] * p1.position[1] * p3.position[2] + p1.position[0] * p2.position[1] * p3.position[2]);
    }
    
    return abs(volume) / 1000;
}

/*
* Calculates surface area of the mesh 
*/
double find_surface_area(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices) 
{
    double sa = 0.0; 

    for (int i = 0; i < triangles.size(); i++) {
        glm::uvec3 triangle = triangles[i];

        glm::vec3 a = vertices[triangle[1]].position - vertices[triangle[0]].position; 
        glm::vec3 b = vertices[triangle[2]].position - vertices[triangle[0]].position; 

        sa += glm::length(glm::cross(a, b)); 
    }

    return 0.5 * sa; 
}

// Checks whether a triangle has an obtuse angle 
bool is_obtuse(Vertex& currentVertex, Vertex& j, Vertex& p)
{
    double cosine_i = glm::dot(glm::normalize(j.position - currentVertex.position), glm::normalize(p.position - currentVertex.position));
    double cosine_j = glm::dot(glm::normalize(currentVertex.position - j.position), glm::normalize(p.position - j.position));
    double cosine_p = glm::dot(glm::normalize(j.position - p.position), glm::normalize(currentVertex.position - p.position));

    return cosine_i < 0 || cosine_j < 0 || cosine_p < 0; 
}

// Computes mixed voronoi region area for a triangle
double mixed_voronoi(Vertex& currentVertex, Vertex& j, Vertex& p, Vertex& q) 
{
    // Find cot(alpha)
    glm::vec3 p_i = currentVertex.position - p.position;
    glm::vec3 p_j = j.position - p.position;
    double alphaWeight = glm::dot(p_i, p_j) / glm::length(glm::cross(p_i, p_j));

    // Find cot(beta)
    glm::vec3 q_i = currentVertex.position - q.position;
    glm::vec3 q_j = j.position - q.position;
    double betaWeight = glm::dot(q_i, q_j) / glm::length(glm::cross(q_i, q_j));
    double norm = glm::length(currentVertex.position - j.position);

    return 0.125 * (alphaWeight + betaWeight) * (norm * norm); 
}

/*
* Finds Voronoi Area
*/
double find_voronoi_area(Vertex& currentVertex,
    std::vector<Vertex>& vertices) 
{
    int ringSize = currentVertex.ring.size();
    double A_i = 0; 
    
    for (int k = 0; k < ringSize; k++) {
        // Get opposite vertex
        Vertex j = vertices[currentVertex.ring[k]];
        // Get next vertex
        Vertex p = vertices[currentVertex.ring[(k + 1) % ringSize]];
        // Get previous vertex
        Vertex q = k == 0 ? vertices[currentVertex.ring[ringSize - 1]] : vertices[currentVertex.ring[(k - 1)]];

        // Non-obtuse, can use Voronoi
        if (!is_obtuse(currentVertex, j, p)) {
            A_i += mixed_voronoi(currentVertex, j, p, q); 
        } else { // Don't use Voronoi with obtuse angles
            double cosine_i = glm::dot(glm::normalize(j.position - currentVertex.position), glm::normalize(p.position - currentVertex.position));
            double area_T = 0.5 * glm::length(glm::cross(j.position - currentVertex.position, p.position - currentVertex.position)); 
            
            if (cosine_i < 0) { // Angle at vertex i is obtuse 
                A_i += 0.5 * area_T; 
            } else {
                A_i += 0.25 * area_T; 
            }
        }
    }
    return  A_i; 
}

/*
* Find the Gaussian curvature k_g at a vertex v. 
*/
float find_gaussian_curvature(Vertex& currentVertex, 
    std::vector<Vertex>& vertices, 
    float A_i)
{
    int ringSize = currentVertex.ring.size();
    double sumTheta = 0; 

    for (int k = 0; k < ringSize; k++) {
        // Get opposite vertex
        Vertex j = vertices[currentVertex.ring[k]];
        // Get next vertex
        Vertex p = vertices[currentVertex.ring[(k + 1) % ringSize]];
        
        // Find theta_j
        glm::vec3 ij = j.position - currentVertex.position; 
        glm::vec3 ip = p.position - currentVertex.position; 
       
        // Find radians and add the theta value
        double input = glm::dot(ij, ip) / (glm::length(ij) * glm::length(ip)); 
        sumTheta += glm::acos(input);
    }

    return (2 * M_PI - sumTheta) / A_i; 
}

/*
* This finds the mean curvature at a vertex v - H
 */
glm::vec3 find_mean_curvature(Vertex& currentVertex, 
    std::vector<Vertex>& vertices,
    float A_i) 
{
    int ringSize = currentVertex.ring.size();
    glm::vec3 laPlace = glm::vec3(0.0); 

    for (int k = 0; k < ringSize; k++) {
        // Get opposite vertex 
        Vertex j = vertices[currentVertex.ring[k]]; 
        // Get next vertex 
        Vertex p = vertices[currentVertex.ring[(k + 1) % ringSize]];
        // Get previous vertex 
        Vertex q = k == 0 ? vertices[currentVertex.ring[ringSize - 1]] : vertices[currentVertex.ring[(k - 1)]]; 

        // Find cot(alpha)
        glm::vec3 p_i = currentVertex.position - p.position; 
        glm::vec3 p_j = j.position - p.position; 
        double alphaWeight = glm::dot(p_i, p_j) / glm::length(glm::cross(p_i,p_j)); 

        // Find cot(beta)
        glm::vec3 q_i = currentVertex.position - q.position;
        glm::vec3 q_j = j.position - q.position;
        double betaWeight = glm::dot(q_i, q_j) / glm::length(glm::cross(q_i, q_j)); 

        laPlace += (float)(alphaWeight + betaWeight) * (currentVertex.position - j.position);
    }

    return laPlace; 
}



// https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
// Calculate the normal vector at a triangle
glm::vec3 calculate_surface_normal(glm::uvec3& triangle, std::vector<Vertex>& vertices)
{
    glm::vec3 normal;

    glm::vec3& p1 = vertices[triangle[0]].position;
    glm::vec3& p2 = vertices[triangle[1]].position;
    glm::vec3& p3 = vertices[triangle[2]].position;

    glm::vec3 u = p2 - p1;
    glm::vec3 v = p3 - p1;

    normal.x = u.y * v.z - u.z * v.y;
    normal.y = u.z * v.x - u.x * v.z;
    normal.z = u.x * v.y - u.y * v.x;

    return normal;
}

// Find list of normals of each triangle
std::vector<Ray> find_normals(std::vector<Vertex>& vertices, std::vector<glm::uvec3>& triangles)
{
    // Find the center of the mesh
    std::vector<glm::vec3> positions;
    std::transform(std::begin(vertices), std::end(vertices),
        std::back_inserter(positions),
        [](const Vertex& v) { return v.position; });

    const glm::vec3 center = std::accumulate(std::begin(positions), std::end(positions), glm::vec3(0.0f)) / static_cast<float>(positions.size());

    // Find normals
    std::vector<Ray> normals = {};

    for (glm::uvec3& tri : triangles) {
        Ray r = {};
        glm::vec3 pos = (1.0f / 3) * (vertices[tri[0]].position + vertices[tri[1]].position + vertices[tri[2]].position);

        glm::vec3 I = pos - center;
        glm::vec3 N = calculate_surface_normal(tri, vertices);

        float dot_product = glm::dot(I, N);

        // Not outside of the surface
        if (dot_product < 0) {
            N = -1.0f * N;
        }

        normals.push_back(Ray { pos, glm::normalize(N), 0.1f });
    }
    return normals;
}

glm::vec3 vertex_normal(int i, std::map<int, std::vector<int>>& vertexToTri, std::vector<Ray>& normals) 
{
    glm::vec3 N; 

    // Sum the normals of the surrounding triangles
    for (int t : vertexToTri[i]) {
        N += normals[t].direction; 
    }

    return glm::normalize(N); 
}

/*
* Calculates the global curvature of a mesh
* Sources:
*   1) http://multires.caltech.edu/pubs/diffGeoOps.pdf
    2) http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation
*   3) http://www.geometry.caltech.edu/pubs/DMSB_III.pdf
*/
double find_curvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices, 
    std::map<int, std::vector<int>>& vertexToTri) 
{
    int vertices_count = 0; 
    double curvature = 0.0; 
    double volume = find_volume(triangles, vertices); 
    double k_reg = std::cbrt((4 * M_PI) / (3 * volume));     
    std::vector<Ray> normals = find_normals(vertices, triangles); 

    for (int i = 0; i < vertices.size(); i++) {
        // If the vertex should be excluded
        if (vertices[i].exclude)
            continue; 
        // Retrieve current vertex and voronoi area 
        Vertex& currentVertex = vertices[i]; 
        glm::vec3 n = vertex_normal(i, vertexToTri, normals); 

        // Find Gaussian curvature K_g and mean curvature H 
        double A_i = find_voronoi_area(currentVertex, vertices); 
        double K_g = find_gaussian_curvature(currentVertex, vertices, A_i); 

        // Mean Curvature 
        glm::vec3 laPlace = find_mean_curvature(currentVertex, vertices, A_i);
       // double H = glm::dot(laPlace, n) >= 0 ? (1.0 / (4 * A_i)) * glm::length(laPlace) : (-1.0 / (4 * A_i)) * glm::length(laPlace); 
        double H =  (1.0 / (4 * A_i)) * glm::length(laPlace); 
        std::cout << "H: " << H << std::endl; 
        
        // Correct for round-off errors 
        K_g = std::abs(K_g) < EPS ? 0.0 : K_g; 
        H = std::abs(H) < EPS ? 0.0 : H; 
           
        // Calculate principle curvatures k_1 and k_2 (+ apply thresholding)
        double k1 = H + sqrt(std::max(0.0, H * H - K_g)); 
        double k2 = H - sqrt(std::max(0.0, H * H - K_g));

        // Print for debugging 
        /* std::cout << "Voronoi Area: " << A_i << std::endl; 
        std::cout << "K_g: " << K_g << std::endl;
        std::cout << "H: " << H << std::endl; 
        std::cout << "K1: " << k1 << std::endl;
        std::cout << "K2: " << k2 << std::endl;*/ 
       
        double vertexCurvature = (0.5 * (k1 + k2));   

        currentVertex.setCurvature(vertexCurvature); 
        curvature += vertexCurvature; 
        vertices_count++; 
    }

    // Average the curvature 
    curvature /= vertices_count; 
    
    return curvature; 
}

// Calculates the heat color at a particular vertex using the curvature value
glm::vec3 heat_color_calculation(const Vertex& vertex, double min, double max) {
    // Uniform curvature
    if (max - min < 1e-6)
        return glm::vec3(0.0f, 1.0f, 0.0f); 
    if (vertex.exclude) { // Excluded in curvature calculation 
        return glm::vec3(0.4, 0.4, 0.4);
    }

    glm::vec3 c = glm::vec3(0.0f); 
    float scaledCurvature = (vertex.curvature - min) / (max - min); 

    if (scaledCurvature < 0.25) {
        scaledCurvature *= 4.0;
        c = glm::vec3(0.0, scaledCurvature, 1.0f);
    } else if (scaledCurvature < 0.50) {
        scaledCurvature = (scaledCurvature - 0.25f) * 4.0f;
        c = glm::vec3(0.0f, 1.0f, 1.0f - scaledCurvature);
    } else if (scaledCurvature < 0.75) {
        scaledCurvature = (scaledCurvature - 0.5f) * 4.0f;
        c = glm::vec3(scaledCurvature, 1.0f, 0.0f);
    } else {
        scaledCurvature = (scaledCurvature - 0.75f) * 4.0f;
        c = glm::vec3(1.0f, 1.0f - scaledCurvature, 0.0f);
    }

    return c; 
}

// Find minimum and maximum curvature 
std::pair<double, double> find_min_max(std::vector<Vertex>& vertices) {
    double min = 100000; 
    double max = -100000; 

    for (auto v : vertices) {
        min = std::min(v.curvature, min); 
        max = std::max(v.curvature, max); 
    }

    std::pair<double, double> minMax = {};
    minMax.first = min; 
    minMax.second = max; 

    return minMax;
}

// Determines the heatcolors for the mesh
std::vector<glm::vec3> heat_color(std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices,
    std::map<int, std::vector<int>>& vertexToTri)
{
    std::vector<glm::vec3> colors = {};
    std::pair<double, double> minMax = find_min_max(vertices);

    for (auto v : vertices) {
        colors.push_back(heat_color_calculation(v, minMax.first, minMax.second));
    }

    std::cout << "min: " << minMax.first << std::endl; 
    std::cout << "max: " << minMax.second << std::endl; 
    return colors; 
}

// Scale the mesh vertices to fit into window
void scale_mesh(std::vector<Vertex>& vertices) {
    for (Vertex& v : vertices) {
        v.position = 0.03f * v.position; 
    }
}

// Find the values that correspond to regional curvature 
std::vector<double> find_regional_curvature(std::vector<Vertex>& vertices) 
{
    std::pair<double, int> cn_1 = { 0, 0 }; 
    std::pair<double, int> cn_2 = { 0, 0 }; 
    std::pair<double, int> cn_3 = { 0, 0 }; 
    std::pair<double, int> cn_4 = { 0, 0 }; 
    std::pair<double, int> cn_5 = { 0, 0 }; 
    std::pair<double, int> cn_6 = { 0, 0 }; 

    for (Vertex& v : vertices) {
        // Do not include in finding regional curvature 
        if (v.exclude)
            continue; 

        switch (v.region) {
        case 1:
            cn_1.first += v.curvature; 
            cn_1.second++; 
            continue;
        case 2:
            cn_2.first += v.curvature;
            cn_2.second++; 
            continue;
        case 3:
            cn_3.first += v.curvature;
            cn_3.second++; 
            continue;
        case 4:
            cn_4.first += v.curvature;
            cn_4.second++; 
            continue;
        case 5:
            cn_5.first += v.curvature;
            cn_5.second++; 
            continue;
        case 6:
            cn_6.first += v.curvature;
            cn_6.second++; 
            continue;
        }
    }

    return { cn_1.first / cn_1.second, 
        cn_2.first / cn_2.second, 
        cn_3.first / cn_3.second, 
        cn_4.first / cn_4.second, 
        cn_5.first / cn_5.second,
        cn_6.first / cn_6.second,
    }; 

}

// Center mesh
void center_mesh(std::vector<Vertex>& vertices) {
    std::vector<glm::vec3> positions;
    std::transform(std::begin(vertices), std::end(vertices),
        std::back_inserter(positions),
        [](const Vertex& v) { return v.position; });

    const glm::vec3 center = std::accumulate(std::begin(positions), std::end(positions), glm::vec3(0.0f)) / static_cast<float>(positions.size());

    for (Vertex& v : vertices) {
        v.position -= center; 
    }
}