#include <framework/mesh.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <math.h>
constexpr auto EPS = 1e-8;
constexpr auto M_PI = 3.14159265358979323846; 

/*
 * Used the source: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 * to find the volume of a 3D mesh.
 */
float findVolume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices)
{
    float volume = 0.0f;

    for (int i = 0; i < triangles.size(); i++) {
        glm::uvec3 triangle = triangles[i];

        Vertex& p1 = vertices[triangle[0]];
        Vertex& p2 = vertices[triangle[1]];
        Vertex& p3 = vertices[triangle[2]];

        volume += (1.0 / 6.0) * (-p3.position[0] * p2.position[1] * p1.position[2] + p2.position[0] * p3.position[1] * p1.position[2] + p3.position[0] * p1.position[1] * p2.position[2] - p1.position[0] * p3.position[1] * p2.position[2] - p2.position[0] * p1.position[1] * p3.position[2] + p1.position[0] * p2.position[1] * p3.position[2]);
    }

    return abs(volume) / 1000.0;
}


float findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices) 
{
    float sa = 0.0f; 

    for (int i = 0; i < triangles.size(); i++) {
        glm::uvec3 triangle = triangles[i];

        glm::vec3 a = vertices[triangle[1]].position - vertices[triangle[0]].position; 
        glm::vec3 b = vertices[triangle[2]].position - vertices[triangle[0]].position; 

        sa += glm::length(glm::cross(a, b)); 
    }

    return 0.5f * sa; 
}

/*
* Returns the unit normal vector to vector p. 
*/
glm::vec3 findNormal(glm::vec3& p, glm::vec3& q, glm::vec3& r) 
{
    glm::vec3 normalVector = glm::cross(q - p, r - p); 

    return glm::normalize(normalVector); 
}

float curvatureTriangle(glm::vec3& p1, glm::vec3& p2, glm::vec3& p3)
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

bool alreadyContains(std::vector<int>& neighbours, int v_j)
{
    for (auto i : neighbours) {
        if (v_j == i)
            return true; 
    }
    return false;
}


float findVoronoiArea(Vertex& currentVertex,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) 
{
    float A_i = 0.0f; 
    glm::vec3 refPoint = currentVertex.position; 
    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to.
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        glm::vec3 a = vertices[triangle[1]].position - vertices[triangle[0]].position;
        glm::vec3 b = vertices[triangle[2]].position - vertices[triangle[0]].position;

        A_i += (0.5 * glm::length(glm::cross(a, b))); 
    }

    return (1.0f / 3.0f) * A_i; 
}


/*
* Find the Gaussian curvature k_g at a vertex v. 
*/
float findGaussianCurvature(Vertex& v, std::vector<int>& adjacentPoints,
    std::map<int, std::vector<int>>& vertrexToTri, std::vector<Vertex>& vertices, 
    float A_i)
{
    glm::vec3 refPoint = v.position; 
    float sumAngle = 0.0f; 

    for (int i = 0; i < adjacentPoints.size(); i++) {
        glm::vec3 p = vertices[adjacentPoints[i]].position; 
        float angle = glm::angle(glm::normalize(refPoint), glm::normalize(p));
        sumAngle += angle;
    }

    float num = (2.0f * M_PI) - sumAngle; 

    return num / A_i; 
}


/*
* Given a vertex, returns all points connected to it. 
*/
std::vector<int> findAdjacentVertices(Vertex& currentVertex,
    int v_i,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) 
{
    
    auto refPoint = currentVertex.position; 
    std::vector<glm::vec3> adjacentPoints = {}; 
    std::vector<int> adjacentIndices = {}; 

    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to. 
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        // Only add the vertex, if it isn't equal to current vertex and we didn't add it yet. 
        if (triangle[0] != v_i && !alreadyContains(adjacentIndices, triangle[0]))
            adjacentIndices.push_back(triangle[0]); 
        if (triangle[1] != v_i && !alreadyContains(adjacentIndices, triangle[1]))
            adjacentIndices.push_back(triangle[1]); 
        if (triangle[2] != v_i && !alreadyContains(adjacentIndices, triangle[2]))
             adjacentIndices.push_back(triangle[2]); 
    }

    return adjacentIndices; 
}

/*
 * Find the _two_ triangles vertex i and j are connected to
 */
std::pair<int, int> findSharedTriangles(int v_i, int v_j, std::map<int, std::vector<int>>& vertexToTri)
{
    std::pair<int, int> pair;
    auto tri_i = vertexToTri[v_i];
    auto tri_j = vertexToTri[v_j];
    std::vector<int> matchedTriangles = {};

    for (int i = 0; i < tri_i.size(); i++) {
        int triangle_i = tri_i[i];
        for (int j = 0; j < tri_j.size(); j++) {
            int triangle_j = tri_j[j];

            // If found one triangle, move to next one
            if (triangle_i == triangle_j) {
                matchedTriangles.push_back(triangle_i);
                break;
            }
        }
        // If two triangles found, move on
        if (matchedTriangles.size() > 1)
            break;
    }

    pair.first = matchedTriangles[0];
    pair.second = matchedTriangles[1];

    return pair;
}


float calculateCotangentWeights(int v_i, int v_j, 
    std::map<int, std::vector<int>>& vertexToTri,
    std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices)
{
    // Find two triangles both vertices are connected to 
    auto twoTriangles = findSharedTriangles(v_i, v_j, vertexToTri); 
    glm::uvec3 triangle1 = triangles[twoTriangles.first];
    glm::uvec3 triangle2 = triangles[twoTriangles.second]; 

    // Find other two vertices 
    int v_q;
    int v_k; 

    for (int z = 0; z < 3; z++) {
        if (triangle1[z] != v_i && triangle1[z] != v_j)
            v_q = triangle1[z]; 
        if (triangle2[z] != v_i && triangle2[z] != v_j)
            v_k = triangle2[z]; 
    }

    /* if (v_i == 3) {
        std::cout << "v_i: " << v_i << std::endl; 
        std::cout << "v_j: " << v_j << std::endl; 
        std::cout << "v_q: " << v_q << std::endl; 
        std::cout << "v_k: " << v_k << std::endl; 
    }*/ 

    // Length of edge (v_i - v_j), stays the same
    float a = glm::distance(vertices[v_i].position, vertices[v_j].position); 

    // Calculate alpha_ij 
    float b = glm::distance(vertices[v_i].position, vertices[v_q].position); 
    float c = glm::distance(vertices[v_j].position, vertices[v_q].position);
    float area = 0.5f * glm::length(glm::cross(vertices[v_i].position - vertices[v_q].position, vertices[v_j].position - vertices[v_q].position)); 
    float alpha_ij = (b*b + c*c - a*a) / (4 * area); 

    // Calculate beta_ij 
    float b2 = glm::distance(vertices[v_i].position, vertices[v_k].position);
    float c2 = glm::distance(vertices[v_j].position, vertices[v_k].position);
    float area2 = 0.5f * glm::length(glm::cross(vertices[v_i].position - vertices[v_k].position, vertices[v_j].position - vertices[v_k].position));
    float beta_ij = (b2 * b2 + c2 * c2 - a * a) / (4 * area2);

    return alpha_ij + beta_ij; 
}


/*
* This finds the mean curvature at a vertex v - H
 */
float findMeanCurvature(Vertex& currentVertex, int vIndex, 
    std::vector<int>& adjacentVertices,
    std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices,
    std::map<int, std::vector<int>>& vertexToTri, 
    float A_i) 
{

    glm::vec3 laplaceP = glm::vec3(0.0f);  
    float meanCurvature = 0.0f; 
    
    // Find the triangles this vertex is connected to 
    auto triangleIndices = vertexToTri[vIndex]; 

    for (int i = 0; i < adjacentVertices.size(); i++) {
        // Retrieve opposite vertex
        int v_j = adjacentVertices[i]; 
        Vertex vertex_j = vertices[v_j]; 

        // (cot(alpha_ij) + cot(beta_ij)) * (v_j - v_i) 
        float cotangentWeights = calculateCotangentWeights(vIndex, v_j, vertexToTri, triangles, vertices); 
        laplaceP = laplaceP + cotangentWeights * (vertex_j.position - currentVertex.position); 
    }

    // Scale with Voronoi area  
    laplaceP = (1 / (2.0f * A_i)) * laplaceP; 
    
    // Find the norm and divide it by 2
    meanCurvature = glm::length(laplaceP) / 2.0f;

    return meanCurvature; 
}

float findCurvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri) 
{
    float curvature = 0.0f; 
    for (int i = 0; i < vertices.size(); i++) {
        // Retrieve current vertex, adjancent vertices, and voronoi area 
        auto currentVertex = vertices[i]; 
        std::vector<int> adjacentPoints = findAdjacentVertices(currentVertex, i, vertexToTri[i], triangles, vertices); 
        float A_i = findVoronoiArea(currentVertex, vertexToTri[i], triangles, vertices); 
      
        // Find Gaussian curvature K_g and mean curvature H 
        float K_g = findGaussianCurvature(currentVertex, adjacentPoints, vertexToTri, vertices, A_i); 
        float H = findMeanCurvature(currentVertex, i, adjacentPoints, triangles, vertices, vertexToTri, A_i); 
       
        // Calculate principle curvatures k_1 and k_2 
        float k1 = H + sqrt(H * H - K_g); 
        float k2 = H - sqrt(H * H - K_g); 

        // Debugging prints
        std::cout << "Voronoi Area: " << A_i << std::endl; 
        std::cout << "K_g: " << K_g << std::endl;
        std::cout << "H: " << H << std::endl; 
        std::cout << "K1: " << k1 << std::endl;
        std::cout << "K2: " << k2 << std::endl; 

        curvature += (0.5 * (k1 + k2)); 
    }

    // Average the curvature 
    curvature /= vertices.size(); 
    
    return curvature; 
}
