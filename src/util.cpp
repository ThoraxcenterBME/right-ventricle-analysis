#include <framework/mesh.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/reciprocal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <math.h>
constexpr auto EPS = 1e-6;
constexpr auto M_PI = 3.14159265358979323846; 

/*
 * Used the source: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 * to find the volume of a 3D mesh.
 */
double findVolume(std::vector<glm::uvec3>& triangles, 
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
double findSurfaceArea(std::vector<glm::uvec3>& triangles, 
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

/*
* Returns the unit normal vector to vector p. 
*/
glm::vec3 findNormal(glm::vec3& p, 
    glm::vec3& q, 
    glm::vec3& r) 
{
    glm::vec3 normalVector = glm::cross(q - p, r - p); 

    return glm::normalize(normalVector); 
}

// Deprecated: Currently not in use 
float curvatureTriangle(glm::vec3& p1, 
    glm::vec3& p2, 
    glm::vec3& p3)
{
    glm::vec3 n1 = findNormal(p1, p2, p3); 
    glm::vec3 n2 = findNormal(p2, p1, p3); 
    glm::vec3 n3 = findNormal(p3, p1, p2); 

    // Edge between p1 and p2 
    float c1 = glm::dot(n2 - n1, glm::normalize(p2 - p1)) / glm::length(p2 - p1);
    
    // Between p1 and p3 
    float c2 = glm::dot(n3 - n1, glm::normalize(p3 - p1)) / glm::length(p3 - p1); 

    // Between p2 and p3 
    float c3 = glm::dot(n3 - n2, glm::normalize(p3 - p2)) / glm::length(p3 - p2);

    return (1.0f / 3.0f) * (c1 + c2 + c3); 
}

/*
* Finds Voronoi Area, simple 1/3 of the area of the surrounding neighbouring triangles
*/
double findVoronoiArea(Vertex& currentVertex,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) 
{
    double A_i = 0.0; 
    glm::vec3 refPoint = currentVertex.position; 
    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to.
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        glm::vec3 a = vertices[triangle[1]].position - vertices[triangle[0]].position;
        glm::vec3 b = vertices[triangle[2]].position - vertices[triangle[0]].position;

        A_i += (0.5 * glm::length(glm::cross(a, b))); 
    }

    return (1.0 / 3.0) * A_i; 
}

/*
* Find the Gaussian curvature k_g at a vertex v. 
*/
float findGaussianCurvature(Vertex& currentVertex, 
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
        glm::vec3 ji = j.position - currentVertex.position; 
        glm::vec3 pi = p.position - currentVertex.position; 
       
        // Find radians and add the theta value
        double input = glm::dot(ji, pi) / (glm::length(ji) * glm::length(pi)); 
        sumTheta += glm::acos(input);
    }

    return (2 * M_PI - sumTheta) / A_i; 
}

/*
* This finds the mean curvature at a vertex v - H
 */
double findMeanCurvature(Vertex& currentVertex, 
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
        glm::vec3 ip = currentVertex.position - p.position; 
        glm::vec3 jp = j.position - p.position; 
        double alpha = glm::dot(ip, jp) / (glm::length(ip) * glm::length(jp)); 

        // Find cot(beta)
        glm::vec3 iq = currentVertex.position - q.position;
        glm::vec3 jq = j.position - q.position; 
        double beta = glm::dot(iq, jq) / (glm::length(iq) * glm::length(jq)); 

        laPlace += (float)(glm::cot(alpha) + glm::cot(beta)) * (j.position - currentVertex.position);

    }

    return (1.0 / (4 * A_i)) * glm::length(laPlace); 
}

/*
* Calculates the global curvature of a mesh
* Sources:
*   1) http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation
*   2) http://www.geometry.caltech.edu/pubs/DMSB_III.pdf
*/
double findCurvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices, 
    std::map<int, std::vector<int>>& vertexToTri) 
{
    double curvature = 0.0; 
    for (int i = 0; i < vertices.size(); i++) {
        // Retrieve current vertex and voronoi area 
        Vertex& currentVertex = vertices[i]; 
        double A_i = findVoronoiArea(currentVertex, vertexToTri[i], triangles, vertices); 
      
        // Find Gaussian curvature K_g and mean curvature H 
        double K_g = findGaussianCurvature(currentVertex, vertices, A_i); 
        double H = findMeanCurvature(currentVertex, vertices, A_i); 
     
        // Correct for round-off errors 
        K_g = std::abs(K_g) < EPS ? 0.0 : K_g; 
        H = std::abs(H) < EPS ? 0.0 : H; 
       
        // Calculate principle curvatures k_1 and k_2 
        double k1 = H + sqrt(H * H - K_g); 
        double k2 = H - sqrt(H * H - K_g);

        // Print for debugging 
        std::cout << "Number of triangles: " << currentVertex.ring.size() << std::endl; 
        std::cout << "Voronoi Area: " << A_i << std::endl; 
        std::cout << "K_g: " << K_g << std::endl;
        std::cout << "H: " << H << std::endl; 
        std::cout << "K1: " << k1 << std::endl;
        std::cout << "K2: " << k2 << std::endl;

        double vertexCurvature = (0.5 * (k1 + k2));  

        std::cout << "Curvature " << i << ": " << vertexCurvature << std::endl; 
        currentVertex.setCurvature(vertexCurvature); 
        curvature += vertexCurvature; 
    }

    // Average the curvature 
    curvature /= vertices.size(); 
    
    return curvature; 
}

// Calculates the heat color at a particular vertex using the curvature value
glm::vec3 heatColorCalculation(const Vertex& vertex, double min, double max) {
    // Uniform curvature
    if (max - min < 1e-6)
        return glm::vec3(0.0f, 1.0f, 1.0f); 

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

std::pair<double, double> findMinMax(std::vector<Vertex>& vertices) {
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
std::vector<glm::vec3> heatColor(std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices,
    std::map<int, std::vector<int>>& vertexToTri)
{
    std::vector<glm::vec3> colors = {};
    std::pair<double, double> minMax = findMinMax(vertices);

    for (auto v : vertices) {
        colors.push_back(heatColorCalculation(v, minMax.first, minMax.second));
    }

    return colors; 
}

void scale(std::vector<Vertex>& vertices) {
    for (Vertex& v : vertices) {
        v.position = 0.01f * v.position; 
    }
}