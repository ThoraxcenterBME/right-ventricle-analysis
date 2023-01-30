#include <framework/mesh.h>
#include <iostream>

struct ColorRegion {
    int key; 
    glm::vec3 c1 = glm::vec3(0.7, 0.1, 0.8);
    glm::vec3 c2 = glm::vec3(0.9, 0.9, 0.1);
    glm::vec3 c3 = glm::vec3(0.7, 0.7, 0.8);
    glm::vec3 c4 = glm::vec3(0.9, 0.5, 0.1);
    glm::vec3 c5 = glm::vec3(1.0, 1.0, 1.0);
    glm::vec3 c6 = glm::vec3(1.0, 0.0, 0.0); 
};


double find_volume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
double find_curvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
double find_surface_area(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 
std::vector<glm::vec3> heat_color(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
void scale(std::vector<Vertex>& vertices); 
std::vector<double> find_regional_curvature(std::vector<Vertex>& vertices); 