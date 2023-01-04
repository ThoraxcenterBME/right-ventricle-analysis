#include <framework/mesh.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/reciprocal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <math.h>
constexpr auto EPS = 1e-7;
constexpr auto M_PI = 3.14159265358979323846; 
constexpr auto conversion = 57.295779; 
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
* Finds Voronoi Area, simple 1/3 of the area of the surrounding neighbouring triangles
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
double find_mean_curvature(Vertex& currentVertex, 
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

        // Check if iterations of the ring is correct, only check for one vertex 
        /* if (currentVertex.index == 37) {
            std::cout << "Current vertex j: " << j.index << std::endl; 
            std::cout << "Forward vertex j + 1: " << p.index << std::endl; 
            std::cout << "Backward vertex j - 1: " << q.index << std::endl; 
        }*/ 

        // Find cot(alpha)
        glm::vec3 p_i = currentVertex.position - p.position; 
        glm::vec3 p_j = j.position - p.position; 
        double alpha = glm::acos(glm::dot(p_i, p_j) / (glm::length(p_i) * glm::length(p_j))); 
        double alphaWeight = glm::dot(p_i, p_j) / glm::length(glm::cross(p_i,p_j)); 

        // Find cot(beta)
        glm::vec3 q_i = currentVertex.position - q.position;
        glm::vec3 q_j = j.position - q.position;
        double beta = glm::acos(glm::dot(q_i, q_j) / (glm::length(q_i) * glm::length(q_j))); 
        double betaWeight = glm::dot(q_i, q_j) / glm::length(glm::cross(q_i, q_j)); 

        glm::vec3 check = j.position - currentVertex.position; 
        laPlace += (float)(alphaWeight + betaWeight) * (currentVertex.position - j.position);
      
        // Some information
        /* if (ringSize == 6 && currentVertex.index == 37) {
            std::cout << "angles: " << alpha*conversion << " " << beta*conversion << std::endl; 
            std::cout << "weights (summed): " << glm::cot(alpha) + glm::cot(beta) << std::endl; 
            std::cout << "vector: " << check.x << " " << check.y << " " << check.z << std::endl; 
        } */ 

    }

    return (1.0 / (4 * A_i)) * glm::length(laPlace); 
}

/*
* Calculates the global curvature of a mesh
* Sources:
*   1) http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation
*   2) http://www.geometry.caltech.edu/pubs/DMSB_III.pdf
*/
double find_curvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices, 
    std::map<int, std::vector<int>>& vertexToTri) 
{
    int vertices_count = 0; 
    double curvature = 0.0; 
    for (int i = 0; i < vertices.size(); i++) {
        // If the vertex should be excluded
        if (vertices[i].exclude)
            continue; 
        // Retrieve current vertex and voronoi area 
        Vertex& currentVertex = vertices[i]; 

        // Find Gaussian curvature K_g and mean curvature H 
        double A_i = find_voronoi_area(currentVertex, vertices); 
        double K_g = find_gaussian_curvature(currentVertex, vertices, A_i); 
        double H = find_mean_curvature(currentVertex, vertices, A_i); 
       
     
        // Correct for round-off errors 
        K_g = std::abs(K_g) < EPS ? 0.0 : K_g; 
        H = std::abs(H) < EPS ? 0.0 : H; 
           
        // Calculate principle curvatures k_1 and k_2 
        double k1 = H + sqrt(std::max(0.0, H * H - K_g)); 
        double k2 = H - sqrt(std::max(0.0, H * H - K_g));

        // Print for debugging 
        /* std::cout << "Number of triangles: " << currentVertex.ring.size() << std::endl; 
        std::cout << "Voronoi Area: " << A_i << std::endl; 
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

void scale(std::vector<Vertex>& vertices) {
    for (Vertex& v : vertices) {
        v.position = 0.03f * v.position; 
        v.position.z += 2.0f; 
    }
}