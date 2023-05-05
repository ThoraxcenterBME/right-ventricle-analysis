#include <framework/mesh.h>
#include <framework/ray.h>
#include <iostream>


struct ColorRegion {
    // List of colours for coloring each region
    std::vector<glm::vec3> colors = { 
        glm::vec3(0.7, 0.1, 0.8), 
        glm::vec3(0.9, 0.9, 0.1), 
        glm::vec3(0.7, 0.7, 0.8), 
        glm::vec3(0.9, 0.5, 0.1), 
        glm::vec3(1.0, 1.0, 1.0), 
        glm::vec3(1.0, 0.0, 0.0),
        glm::vec3(0.0, 0.0, 0.0)
    };
};

// Info for displaying on GUI
struct RVInfo {
    float volume;
    float surfaceArea;
    float curvature;
    float radius; 
    std::vector<float> regional_curvs; 
    double global_area_strain;
};

// Info that can be outputed to CSV
struct PrintInfo {
    float volume;
    float surface_area;
    float curvature;
    float index_curv;
    float min_curv;
    float max_curv;
    std::vector<double> curvatures;
    std::vector<double> volumes;
    std::vector<double> surface_areas;
    std::pair<double, double> minmax;
};

double find_volume(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices);
double find_curvature(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
double find_surface_area(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices); 
std::vector<glm::vec3> heat_color(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
void scale_mesh(std::vector<Vertex>& vertices); 
std::vector<double> find_regional_curvature(std::vector<Vertex>& vertices); 
void center_mesh(std::vector<Vertex>& vertices); 
std::vector<Ray> find_normals(std::vector<Vertex>& vertices, std::vector<glm::uvec3>& triangles); 
std::vector<Ray> findLaplaceRays(std::vector<glm::uvec3>& triangles, std::vector<Vertex>& vertices, std::map<int, std::vector<int>>& vertexToTri);
std::vector<double> regional_volumes(std::vector<Vertex>& vs, std::vector<glm::uvec3>& ts, std::map<int, std::vector<int>>& vertexToTri);
std::vector<double> regional_surface_areas(std::vector<Vertex>& vs, std::vector<glm::uvec3>& ts, std::map<int, std::vector<int>>& vertexToTri); 
double find_indexed_curvature(std::vector<Vertex>& vertices);
std::pair<double, double> find_min_max(std::vector<Vertex>& vertices); 
