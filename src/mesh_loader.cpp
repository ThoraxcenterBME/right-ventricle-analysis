#include "mesh_loader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <glm/gtx/normal.hpp>

void setVertexPosition(Vertex& v, float x, float y, float z) {
    glm::vec3 pos = glm::vec3(x, y, z); 
    v.position = pos; 
}

// Determine which vertex is not yet explored
int findExcludedVertex(int i, int j, glm::uvec3& tri) {
    for (int k = 0; k < 3; k++) {
        if (tri[k] != i && tri[k] != j) {
            return tri[k]; 
        }
    }
}

// Find the two triangles shared by two vertices
std::pair<int, int> findSharedTriangle(std::vector<int>& triA, std::vector<int>& triB) {
    std::pair<int, int> triangles; 
    bool foundFirst = false; 

    for (int i = 0; i < triA.size(); i++) { 
        for (int j = 0; j < triB.size(); j++) {
            if (triA[i] == triB[j] && !foundFirst) {
                triangles.first = triA[i];
                foundFirst = true; // set flag
            } else if (triA[i] == triB[j]) {
                triangles.second = triA[i];
            }
        }
    }
    return triangles; 
}

// Returns which triangle to visit next 
// 0 - triangle A  (first)
// 1 - triangle B  (second)
int findNextTriangle(glm::uvec3& triangleA, glm::uvec3& triangleB, std::vector<int>& ring) {
    // Visit the triangle that has least visited vertices 
    int countA = 0; 
    int countB = 0;

    for (int k = 0; k < 3; k++) {
        int vA = triangleA[k]; 
        int vB = triangleB[k]; 

        for (int r = 0; r < ring.size(); r++) {
            if (ring[r] == vA)
                countA++; 
            if (ring[r] == vB)
                countB++; 
        }
    }

    return countA <= countB ? 0 : 1; 
}

int findAdjacentVertex(glm::uvec3& tri, int i) {
    for (int k = 0; k < 3; k++) {
        if (tri[k] != i) {
            return tri[k]; 
        }
    }
}

void initializeRings(std::map<int, std::vector<int>>& vertexMap, 
    std::vector<Vertex>& vertices, 
    std::vector<glm::uvec3>& triangles) {
    
    for (int i = 0; i < vertices.size(); i++) { 
        // Instantiate 1-neighbourhood ring
        std::vector<int> ring = {}; 
        int vertexCount = vertexMap[i].size(); // Number of vertices in ring
        int nextVertex = findAdjacentVertex(triangles[vertexMap[i][0]], i); // Start at an arbitrary vertex , gives you a triangle index, not vertex index

        // Continue to find ring vertices until ring is complete 
        while (ring.size() < vertexCount) {
            ring.push_back(nextVertex); 

            // If ring is complete, stop 
            if (ring.size() == vertexCount)
                break; 

            std::pair<int, int> sharedTri = findSharedTriangle(vertexMap[i], vertexMap[nextVertex]);
            int nextTri = findNextTriangle(triangles[sharedTri.first], triangles[sharedTri.second], ring); 

            // Find the next vertex to visit
            if (nextTri == 0) {
                nextVertex = findExcludedVertex(i, nextVertex, triangles[sharedTri.first]); 
            } else {
                nextVertex = findExcludedVertex(i, nextVertex, triangles[sharedTri.second]);   
            }
         
        }

        // Set the ring of vertices
        vertices[i].ring = ring; 
    }
}

// Function for printing vertex map
void printVertexMap(std::map<int, std::vector<int>>& vmap) 
{
    for (int i = 0; i < vmap.size(); i++) {
        std::cout << "Vertex: " << i << std::endl; 
        for (int j = 0; j < vmap[i].size(); j++) {
            std::cout << "triangle: " << vmap[i][j] << std::endl; 
        }
    }
}

// Function for printing vertex map
void printVertexRing(std::vector<Vertex>& vertices)
{

    for (int i = 0; i < vertices.size(); i++) {
        std::cout << vertices[i].ring.size() << " " << std::flush;
        for (int j = 0; j < vertices[i].ring.size(); j++) {
            std::cout << vertices[i].ring[j] << " " << std::flush;
        }
        std::cout << "" << std::endl;
    }
}

// Load the ring from a file
void loadRingFromFile(std::string& name, std::vector<Vertex>& vertices) {
    std::string line;
    std::ifstream myfile;
    myfile.open(std::filesystem::path(DATA_DIR) / name);
    int vIndex = 0; 
 
    while (std::getline(myfile, line)) {
        std::istringstream lines(line);
        
        int neighbourCount; 
        lines >> neighbourCount;

        for (int i = 0; i < neighbourCount; i++) {
            int nn; 
            lines >> nn; 
            vertices[vIndex].ring.push_back(nn); 
        }
        vIndex++; 
    }

}

double largestDist(std::vector<Vertex>& vertices) {
    float dist = -10000; 
    for (int i = 0; i < vertices.size(); i++) {
        for (int j = 0; j < vertices.size(); j++) {
            dist = std::max(dist, glm::length(vertices[i].position - vertices[j].position)); 
        }
    }
    return dist; 
}

// Mark some vertices as not included in curvature calculation 
void mark_excluded(std::string& name, std::vector<Vertex>& vertices)
{
    std::ifstream myfile;
    myfile.open(std::filesystem::path(DATA_DIR) / name);
    std::string linebuf; 

    while (std::getline(myfile, linebuf)) {
        std::istringstream lines(linebuf);
        int key;
        lines >> key; 
        vertices[key].exclude = true; 
    }
}

void mark_regions(std::string& in, std::vector<Vertex>& vertices)
{   
    std::ifstream myfile;
    myfile.open(std::filesystem::path(DATA_DIR) / in);
    std::string linebuf;

    while (std::getline(myfile, linebuf)) {
        std::istringstream lines(linebuf);
        std::string lineType;
        lines >> lineType;
        int key; 
        lines >> key; 

        // Inflow Tract = 1 
        if (lineType == "it") {
            vertices[key].set_region(1); 
            continue; 
        }
        // Outflow Tract = 2
        if (lineType == "ot") {
            vertices[key].set_region(2);
            continue; 
        }
        // Septal Body = 3
        if (lineType == "sb") {
            vertices[key].set_region(3);
            continue; 
        }
        // Free wall Body = 4 
        if (lineType == "fb") {
            vertices[key].set_region(4);
            continue; 
        }
        // Septal Apex = 5
        if (lineType == "sa") {
            vertices[key].set_region(5);
            continue; 
        }
        // Free-wall Apex = 6 
        if (lineType == "fa") {
            vertices[key].set_region(6);
            continue; 
        }
    }
}

Mesh loadMeshRV(std::istream& in) {

    Mesh rv = {}; 
    std::vector<Vertex> vertices = {}; 
    std::vector<glm::uvec3> triangles = {}; 
    std::map<int, std::vector<int>> vertexMap = {}; 
    int vertexKey = 0; 
    int triangleIndex = 0; 
    std::string linebuf;
    double minCoordinate = 10000; 
    double maxCoordinate = -10000; 
   
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
             v.index = vertexKey; 
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
    rv.radius = largestDist(rv.vertices) / 2.0; 

    // Printing the rings or need them for the first time
    // initializeRings(rv.vertexToTri, rv.vertices, rv.triangles); 
    // printVertexRing(rv.vertices); 
    
    return rv; 
}

