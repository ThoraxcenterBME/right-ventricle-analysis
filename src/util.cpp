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


float findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices) {
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
glm::vec3 findNormal(glm::vec3& p, glm::vec3& q, glm::vec3& r) {
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

bool alreadyContains(std::vector<glm::vec3>& vertices, glm::vec3& point, glm::vec3& refPoint)
{
    for (auto v : vertices) {
        // If adjacent points already contains it or if it's reference point
        if (glm::distance(v, point) < EPS || glm::distance(v, refPoint) < EPS) {
            return true;
        }
    }
    return false;
}

std::vector<glm::vec3> constructRing(Vertex& currentVertex,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) {
    // Declare the ring 
    std::vector<glm::vec3> ring = {}; 


    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to.
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        glm::vec3 p1 = vertices[triangle[0]].position; 
        glm::vec3 p2 = vertices[triangle[1]].position; 
        glm::vec3 p3 = vertices[triangle[2]].position; 


    }
    return ring;

}

float findVoronoiArea(Vertex& currentVertex,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) {

    float A_i = 0.0f; 
    glm::vec3 refPoint = currentVertex.position; 
    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to.
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        glm::vec3 a = vertices[triangle[1]].position - vertices[triangle[0]].position;
        glm::vec3 b = vertices[triangle[2]].position - vertices[triangle[0]].position;

        A_i += glm::length(glm::cross(a, b)); 
    }

    return (1.0f / 3.0f) * A_i; 
}


/*
* Find the Gaussian curvature k_g at a vertex v. 
*/
float findGaussianCurvature(Vertex& v, std::vector<glm::vec3> adjacentPoints, float A_i)
{
    glm::vec3 refPoint = v.position; 
    float sumAngle = 0.0f; 

    for (auto p : adjacentPoints) {
        float angle = glm::angle(glm::normalize(refPoint), glm::normalize(p));
        sumAngle += angle;
    }

    float num = (2.0f * M_PI) - sumAngle; 

    return num / A_i; 
}

int findPrevPoint() {
    return 0; 
}



/*
* Given a vertex, returns all points 
* Note that this also reshuffles the adjacent vertices, so 
* that all vertices form a _ring_ around the reference vertex.
*/
std::vector<glm::vec3> findAdjacentVertices(Vertex& currentVertex,
    std::vector<int>& triangleIndices,
    std::vector<glm::uvec3>& triangles,
    std::vector<Vertex>& vertices) {
    
    auto refPoint = currentVertex.position; 
    std::vector<glm::vec3> adjacentPoints = {}; 

    for (int i = 0; i < triangleIndices.size(); i++) {
        // Find triangle vertex is connected to. 
        glm::uvec3 triangle = triangles[triangleIndices[i]];

        // Find all points of the triangle
        glm::vec3 p1 = vertices[triangle[0]].position;
        glm::vec3 p2 = vertices[triangle[1]].position;
        glm::vec3 p3 = vertices[triangle[2]].position;

        // Add adjacent points that did not exist already
        if (!alreadyContains(adjacentPoints, p1, refPoint))
            adjacentPoints.push_back(p1); 
        if (!alreadyContains(adjacentPoints, p2, refPoint))
            adjacentPoints.push_back(p2); 
        if (!alreadyContains(adjacentPoints, p3, refPoint))
             adjacentPoints.push_back(p3); 
    }



    return adjacentPoints; 
}

/*
* This finds the curvature of the mesh
 */
float findMeanCurvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices) {
    float curvature = 0.0f; 
    
    // Sum the curvature at each vertex 
    for (int i = 0; i < triangles.size(); i++) {
        glm::uvec3 triangle = triangles[i];

        auto p1 = vertices[0].position; 
        auto p2 = vertices[1].position; 
        auto p3 = vertices[2].position; 

        curvature += curvatureTriangle(p1, p2, p3); 

    }

    // average the curvature value over all vertices 
    float numTriangles = (float)triangles.size(); 

    return curvature / numTriangles; 
}
