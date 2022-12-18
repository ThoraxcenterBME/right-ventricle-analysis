#include <framework/mesh.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/reciprocal.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <map>
#include <math.h>
constexpr auto EPS = 1e-6;
constexpr auto M_PI = 3.14159265358979323846; 

// TODO: make ring around current vertex, clockwise order

/*
 * Used the source: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
 * to find the volume of a 3D mesh.
 */
double findVolume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices)
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


double findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices) 
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
glm::vec3 findNormal(glm::vec3& p, glm::vec3& q, glm::vec3& r) 
{
    glm::vec3 normalVector = glm::cross(q - p, r - p); 

    return glm::normalize(normalVector); 
}

// Deprecated: Currently not in use 
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

// Assign the vertices of triangles, so that it isn't the same as the reference point
void assignVertices(glm::uvec3& triangle, const int v_i, glm::vec3& a, 
    glm::vec3& b, std::vector<Vertex>& vertices) 
{
    bool setA = false; 
    bool setB = false; 

    if (triangle[0] != v_i) {
        a = vertices[triangle[0]].position;
        setA = true; 
    }
    if (triangle[1] != v_i && setA) {
        b = vertices[triangle[1]].position;
        setB = true; 
    } else if (triangle[1] != v_i && !setA) {
        a = vertices[triangle[1]].position; 
        setA = true; 
    }
    if (triangle[2] != v_i && !setB) {
        b = vertices[triangle[2]].position;
        setB = true; 
    }
}


/*
* Find the Gaussian curvature k_g at a vertex v. 
*/
float findGaussianCurvature(Vertex& v, int v_i, std::vector<int>& adjacentPoints,
    std::map<int, std::vector<int>>& vertrexToTri, std::vector<Vertex>& vertices, 
    std::vector<glm::uvec3>& triangles, float A_i)
{
    glm::vec3 refPoint = v.position; 
    std::vector<int> connectedTri = vertrexToTri[v_i]; 
    float sumAngle = 0.0f; 

    for (int i = 0; i < connectedTri.size(); i++) {
        // The other two vertices of the same triangle 
        glm::vec3 a; 
        glm::vec3 b; 
        glm::uvec3 tri = triangles[connectedTri[i]]; 

        assignVertices(tri, v_i, a, b, vertices); 

        glm::vec3 p = a - refPoint; 
        glm::vec3 q = b - refPoint; 

        // To prevent divisions by zero
        if (glm::length(p) < EPS || glm::length(q) < EPS) {
            std::cout << "Division by zero encountered in k_G" << std::endl; 
            return 0.0; 
        }

        double input = glm::dot(p, q) / (glm::length(p) * glm::length(q)); 

        double angle = glm::acos(input); 
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

double findCurvature(std::vector<glm::uvec3>& triangles, 
    std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri) 
{
    double curvature = 0.0; 
    for (int i = 0; i < vertices.size(); i++) {
        // Retrieve current vertex, adjancent vertices, and voronoi area 
        auto currentVertex = vertices[i]; 
        std::vector<int> adjacentPoints = findAdjacentVertices(currentVertex, i, vertexToTri[i], triangles, vertices);  
        double A_i = findVoronoiArea(currentVertex, vertexToTri[i], triangles, vertices); 
      
        // Find Gaussian curvature K_g and mean curvature H 
        double K_g = findGaussianCurvature(currentVertex, i, adjacentPoints, vertexToTri, vertices, triangles, A_i); 
        double H = findMeanCurvature(currentVertex, vertices, A_i); 
     
        // Correct for round-off errors 
        K_g = K_g < EPS && K_g > 0 ? 0.0 : K_g; 
        H = H < EPS && H > 0 ? 0.0 : H; 
       
        // Calculate principle curvatures k_1 and k_2 
        double k1 = H + sqrt(H * H - K_g); 
        double k2 = H - sqrt(H * H - K_g);

        std::cout << "Number of triangles: " << adjacentPoints.size() << std::endl; 
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
