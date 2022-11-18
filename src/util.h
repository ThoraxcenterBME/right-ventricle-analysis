#include <framework/mesh.h>
#include <iostream>

float findVolume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
float findCurvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
float findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 