#include <framework/mesh.h>
#include <iostream>

float findVolume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
float findMeanCurvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
float findSurfaceArea(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 