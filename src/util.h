#include <framework/mesh.h>
#include <iostream>

double findVolume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
double findCurvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
double findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 