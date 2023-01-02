#include <framework/mesh.h>
#include <iostream>

double find_volume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
double find_curvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
double find_surface_area(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 
std::vector<glm::vec3> heat_color(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
void scale(std::vector<Vertex>& vertices); 