#include "mesh_loader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <glm/gtx/normal.hpp>
//struct Mesh {
//    // Vertices contain the vertex positions and normals of the mesh.
//    std::vector<Vertex> vertices;
//    // A triangle contains a triplet of values corresponding to the indices of the 3 vertices in the vertices array.
//    std::vector<glm::uvec3> triangles;
//
//    Material material;
//};

//struct Vertex {
//    glm::vec3 position;
//    glm::vec3 normal;
//    glm::vec2 texCoord; // Texture coordinate
//
//    [[nodiscard]] constexpr bool operator==(const Vertex&) const noexcept = default;
//};

void setVertexPosition(Vertex& v, float x, float y, float z) {
    glm::vec3 pos = glm::vec3(x, y, z); 
    v.position = pos; 
}

void addVertexTriangleMap(std::map<int, std::vector<int>>& vertexMap, 
    std::vector<Vertex>& vertices, 
    std::vector<glm::uvec3>& triangles) {

}

Mesh loadMeshRV(std::istream& in) {
    Mesh rv = {}; 
    std::vector<Vertex> vertices = {}; 
    std::vector<glm::uvec3> triangles = {}; 
    std::map<int, std::vector<int>> vertexMap = {}; 
    int vertexKey = 0; 
    int triangleIndex = 0; 
    std::string linebuf;
   
    while (std::getline(in, linebuf)) { 
        std::istringstream lines(linebuf);
        std::string lineType;
        lines >> lineType;

        // Vertex
        if (lineType == "v") {
            // Load in vertex 
            float x, y, z = 0.0f;
            lines >> x >> y >> z;
            Vertex v = {}; 
            setVertexPosition(v, x, y, z); 

            // Add vertex to list of vertices
            vertices.push_back(v); 
            // Add an empty list entry to the map 
            vertexMap[vertexKey] = {}; 
            // Update vertex key value 
            vertexKey = vertexKey + 1; 
            continue; 
        }
        // Face (Triangle)
        if (lineType == "f") {
            // Load in triangle
            int v1, v2, v3;
            lines >> v1 >> v2 >> v3;
            glm::uvec3 tri = glm::uvec3(v1 - 1, v2 - 1, v3 - 1); 

            // Add triangles to list of triangles 
            triangles.push_back(tri); 
            // Add the triangle to vertex map 
            vertexMap[v1 - 1].push_back(triangleIndex); 
            vertexMap[v2 - 1].push_back(triangleIndex);
            vertexMap[v3 - 1].push_back(triangleIndex); 
            // Update the triangle index 
            triangleIndex = triangleIndex + 1; 
            continue; 
        }
        // If doesn't match anything, continue
        continue; 
    }

    // Set triangles, vertices, and the map 
    rv.triangles = triangles; 
    rv.vertices = vertices; 
    rv.vertexToTri = vertexMap; 

    return rv; 
}

