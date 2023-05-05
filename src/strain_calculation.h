#include <framework/mesh.h>
#include <framework/ray.h>
#include <iostream>

struct Strain {
    // Global value
    double global_es_area;
    double global_ed_area;
    double global_area_strain;

    // List of areas for regions
    std::vector<double> ed_areas;
    std::vector<double> es_areas;
    std::vector<double> strain_values;

    // Three axes for which strain is defined. 
    glm::vec3 long_axis;
    glm::vec3 circ_axis; 
    glm::vec3 radial_axis;
    // Global strain values 
    double longitudinal_strain;
    double circumferential_strain;
    double radial_strain; 

    std::vector<Vertex> vertices_ed;
    std::vector<Vertex> vertices_es;
    std::vector<double> l_strain_values;
    std::vector<double> c_strain_values; 
    std::vector<double> r_strain_values; 
};


double area_strain(double ed_area, double es_area);
void set_regional_area_strain(Strain& strain);
double longitudinal_strain(std::vector<Vertex>& vertices_es,
    std::vector<Vertex>& vertices_ed,
    std::vector<glm::uvec3> triangles,
    std::map<int, std::vector<int>>& vertexToTri,
    glm::vec3& l_axis, Strain& strain);
glm::vec3 find_long_axis(std::vector<Vertex>& vertices, int lower_point, int c1, int c2);