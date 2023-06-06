#include "draw.h"
#include "util.h"
#include "shade.h"
#include "mesh_loader.h"
#include "strain_calculation.h"
// Disable warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/vector_relational.hpp>
#include <imgui/imgui.h>
DISABLE_WARNINGS_POP()
#include <array>
#include <framework/mesh.h>
#include <framework/ray.h>
#include <framework/trackball.h>
#include <framework/window.h>
#include <fstream>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <vector>

bool show_normal_gui = false; 
bool draw_regions_gui = false; 
bool number_vertices = false; 

enum class DiffuseMode {
    None,
    Lambert,
    Gooch
};

enum class SpecularMode {
    None,
    Phong,
    BlinnPhong
};

// Different display modes.
struct DisplayModes {
    bool debug = false;
    DiffuseMode diffuse = DiffuseMode::Lambert;
    SpecularMode specular = SpecularMode::None;
};

enum class LightPlacement {
    Sphere,
    Shadow,
    Specular
};

struct Light {
    glm::vec3 position;
    glm::vec3 color;
};

struct ProgramState {
    DisplayModes displayModes;
    glm::vec3 backgroundColor { 1.0f, 1.0f, 0.878f };

    std::vector<Light> lights;
    int selectedLight = 0;
    LightPlacement lightPlacement = LightPlacement::Sphere;

    Mesh myMesh;
    unsigned selectedVertex = 0xFFFFFFFF;
    bool showSelectedVertex = true;
    MaterialInformation materialInformation;
};

glm::vec3 computeLighting(const ProgramState& programState, unsigned vertexIndex, const glm::vec3& cameraPos, const Light& light)
{
    const auto& vertex = programState.myMesh.vertices[vertexIndex];
    const Positions positions {
        .vertex = vertex.position,
        .light = light.position,
        .camera = cameraPos
    };
    const glm::vec3& normal = vertex.normal;

    // do not change any global variables here! This function will be called for EACH vertex
    // of the mesh, so your change would happen several times
    if (programState.displayModes.debug) {
        return debugColor(programState.materialInformation, positions, normal, light.color);
    }

    glm::vec3 result { 0.0f };
    switch (programState.displayModes.diffuse) {
    case DiffuseMode::Lambert: {
        result += diffuseOnly(programState.materialInformation, positions, normal, light.color);
    } break;
    case DiffuseMode::Gooch: {
        result += gooch(programState.materialInformation, positions, normal, light.color);
    } break;
    case DiffuseMode::None:
        break;
    };

    const glm::vec3 diff = light.position - vertex.position;
    const float dist2 = glm::dot(diff, diff);
    return result / dist2;
}

static std::optional<glm::vec3> getWorldPositionOfPixel(const Window& window, const Trackball& trackball, const glm::vec2& pixel);

static size_t getClosestVertexIndex(const Mesh& mesh, const glm::vec3& pos)
{
    const auto iter = std::min_element(
        std::begin(mesh.vertices), std::end(mesh.vertices),
        [&](const Vertex& lhs, const Vertex& rhs) {
            return glm::length(lhs.position - pos) < glm::length(rhs.position - pos);
        });
    return std::distance(std::begin(mesh.vertices), iter);
}

void keyboard(unsigned char key, ProgramState& state, const Window& window, const Trackball& camera)
{
    switch (key) {
    case ' ': {
        const auto worldPoint = getWorldPositionOfPixel(window, camera, window.getCursorPixel());
        if (worldPoint) {
            state.selectedVertex = unsigned(getClosestVertexIndex(state.myMesh, *worldPoint));
        }
    } break;
    };
}

// Get the 3D world position of the mouse cursor (assuming the depth buffer has been filled).
static std::optional<glm::vec3> getWorldPositionOfPixel(const Window& window, const Trackball& trackball, const glm::vec2& pixel)
{
    float depth;
    glReadPixels(static_cast<int>(pixel.x), static_cast<int>(pixel.y), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

    if (depth < 0.0f || depth >= 1.0f) {
        // This is a work around for a bug in GCC:
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80635
        //
        // This bug will emit a warning about a maybe uninitialized value when writing:
        // return {};
        constexpr std::optional<glm::vec3> tmp;
        return tmp;
    }

    // Coordinates convert from pixel space to OpenGL screen space (range from -1 to +1)
    const glm::vec3 win { pixel, depth };

    // View matrix
    const glm::mat4 view = trackball.viewMatrix();
    const glm::mat4 projection = trackball.projectionMatrix();

    const glm::vec4 viewport { glm::vec2(0), window.getFrameBufferSize() };
    return glm::unProject(win, view, projection, viewport);
}

void draw(const ProgramState& state, const Trackball& camera, std::vector<glm::vec3> vertexColors, std::vector<Vertex>& vertices)
{
    glClearColor(state.backgroundColor.r, state.backgroundColor.g, state.backgroundColor.b, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Initialize projection and view matrices.
    glMatrixMode(GL_PROJECTION);
    const auto projectionMatrix = camera.projectionMatrix();
    glLoadMatrixf(glm::value_ptr(projectionMatrix));
    glMatrixMode(GL_MODELVIEW);
    const auto viewMatrix = camera.viewMatrix();
    glLoadMatrixf(glm::value_ptr(viewMatrix));

    // Drawing mode options.
    glEnable(GL_DEPTH_TEST); // Enable depth test.
    glDepthMask(GL_TRUE); // Enable depth write.
    glDisable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);

    drawMeshWithColors(state.myMesh, vertexColors);
  //  drawAxis();

    glm::vec3 l_axis = find_long_axis(vertices, 906, 102, 63); 
    glm::vec3 r_axis = find_radial_axis(vertices, l_axis); 
    glm::vec3 c_axis = find_circumferential_axis(vertices, l_axis, r_axis); 

    drawRay(Ray { vertices[906].position, l_axis, glm::length(l_axis) }, glm::vec3(1, 0, 0)); 
    drawRay(Ray { glm::vec3(0, 0, 0), r_axis, glm::length(l_axis) }, glm::vec3(0, 1, 0)); 
    drawRay(Ray { glm::vec3(0, 0, 0), c_axis, glm::length(l_axis) }, glm::vec3(0, 0, 1)); 

    // Disable depth write because we don't want the points in our depth buffer (messes with user interaction).
    glDepthMask(GL_FALSE); // Disable depth write.

    // Draw a big yellow point (square) around the selected light.
    glPointSize(40);
    glColor3f(1, 1, 0);
    glBegin(GL_POINTS);
    glVertex3fv(glm::value_ptr(state.lights[static_cast<size_t>(state.selectedLight)].position));
    glEnd();

    // Draw lights as points (squares) in the lights color.
    glPointSize(10);
    glBegin(GL_POINTS);
    for (const auto& light : state.lights) {
        glColor3fv(glm::value_ptr(light.color));
        glVertex3fv(glm::value_ptr(light.position));
    }
    glEnd();

    // Draw a small red point (square) at the selected vertex.
    if (state.showSelectedVertex && state.selectedVertex != 0xFFFFFFFF) {
        glBegin(GL_POINTS);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3fv(glm::value_ptr(state.myMesh.vertices[state.selectedVertex].position));
        glEnd();
    }

    glDepthMask(GL_TRUE); // Disable depth write.
}

void computeLighting(const ProgramState& state, const glm::vec3& cameraPos, std::span<glm::vec3> outVertexColors)
{
    for (unsigned v = 0; v < state.myMesh.vertices.size(); v++) {
        outVertexColors[v] = glm::vec3(0.0f);
        for (const auto& light : state.lights)
            outVertexColors[v] += computeLighting(state, v, cameraPos, light);
    }
}


void drawUI(ProgramState& state, const Trackball& camera, RVInfo& info, 
    std::vector<Ray>& normals)
{
    ColorRegion color = {}; 
    ImGui::Begin("View RV Data");

    // Display Volume
    std::string volumeString = "Volume: " + std::to_string(info.volume);
    ImGui::Text(volumeString.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Surface Area
    std::string saString = "Surface Area: " + std::to_string(info.surfaceArea);
    ImGui::Text(saString.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Radius 
    std::string radString = "Radius: " + std::to_string(info.radius);
    ImGui::Text(radString.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Curvature
    std::string curveString = "Indexed Curvature: " + std::to_string(info.curvature);
    ImGui::Text(curveString.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Check for drawing normals
    ImGui::Checkbox("Show Surface Normals", &show_normal_gui); 
    if (show_normal_gui) {
        auto red = glm::vec3(1, 0, 0); 
        draw_rays(normals, red); 
    }
    ImGui::Spacing();
    ImGui::Separator();

    // Show regions
    ImGui::Checkbox("Color Regions", &draw_regions_gui);
    if (draw_regions_gui) {
        draw_regions(state.myMesh.vertices);
    }
    ImGui::Spacing();
    ImGui::Separator();


    // Display Region 1 
    std::string r1_curvature = "Curvature (Inferior Free Wall): " + std::to_string(info.regional_curvs[0]);
    ImGui::TextColored(ImVec4(color.colors[0].x, color.colors[0].y, color.colors[2].z, 1.0), r1_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 2
    std::string r2_curvature = "Curvature (Lateral Free Wall): " + std::to_string(info.regional_curvs[1]);
    ImGui::TextColored(ImVec4(color.colors[1].x, color.colors[1].y, color.colors[1].z, 1.0), r2_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 3
    std::string r3_curvature = "Curvature (Anterior Free Wall): " + std::to_string(info.regional_curvs[2]);
    ImGui::TextColored(ImVec4(color.colors[2].x, color.colors[2].y, color.colors[2].z, 1.0), r3_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 4
    std::string r4_curvature = "Curvature (Septal Body): " + std::to_string(info.regional_curvs[3]);
    ImGui::TextColored(ImVec4(color.colors[3].x, color.colors[3].y, color.colors[3].z, 1.0), r4_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    ImGui::Text("Red is longitudinal axis");
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Green is radial axis"); 
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Blue is circumferential axis"); 

    ImGui::End();
}

void printHelp()
{
    std::cout << "Program Usage:" << std::endl;
    std::cout << "SPACE - replaces mouse click for selection, will then call your light placement function" << std::endl;
}

// Change 
void set_regional(RVInfo& info, std::vector<Vertex>& vertices) {
    std::vector<double> curvatures = find_regional_curvature(vertices); 

    for (auto& c : curvatures) {
        info.regional_curvs.push_back(c);
    }
}

void write_to_file(std::string filename, PrintInfo info, int i) {
    std::ofstream datafile; 
    datafile.open(std::filesystem::path(DATA_DIR) / filename, std::ios_base::app); 
    datafile << i << ", "; 
    datafile << info.volume << ", ";
    datafile << info.surface_area << ", "; 
    datafile << info.index_curv << ", "; 
   

    for (auto& v : info.volumes) {
        datafile << v / 1000.0 << ", "; 
    }

    for (auto& sa : info.surface_areas) {
        datafile << sa << ", "; 
    }

    for (auto& c : info.curvatures) {
        datafile << c << ", "; 
    }
    datafile << info.minmax.first << ", "; 
    datafile << info.minmax.second; 

    datafile << "\n"; 

    datafile.close(); 
}

std::string construct_file_string(int n)
{
    // healthy-1: Young healthy volunteer_00 (38)
    // healthy-2: RV beutel young healthy volunteer BPD 7_00 (65)
    // healthy-3: RV beutel young healthy volunteer BPD 15_00 (51)
    // healthy-4: RV beutel young healthy volunteer BPD 16_00 (37)
    // healthy-5: RV beutel young healthy volunteer BPD 31_00 (78)
    // tof-1: ToF RV Beutel_00 (41)
    // tof-2: ToF_2_00 (21)
    // asd-1: Post_ASD_1_00 (34)
    if (n / 10 < 1) {
        return "healthy-1/Young healthy volunteer_00" + std::to_string(n) + ".obj";
    } else {
        return "healthy-1/Young healthy volunteer_0" + std::to_string(n) + ".obj";
    }
}

/*
 * Method for postprocessing RV beutel mesh
 */
Mesh get_mesh(std::string& filename)
{
    // Data files that encode regional RV information
    std::string ring = "ring-indices.txt"; // ring-indices, ring-sphere ring-large
    std::string exclude_vertices = "exclude.txt";
    std::string regions = "region-v2.txt";

    std::ifstream ifile;
    ifile.open(std::filesystem::path(DATA_DIR) / filename);
    Mesh rv = loadMeshRV(ifile);
    loadRingFromFile(ring, rv.vertices);

    // Extra post-processing for RV regions
    center_mesh(rv.vertices);
    mark_excluded(exclude_vertices, rv.vertices);
    mark_regions(regions, rv.vertices);

    return rv;
}

void write_strain_to_file(std::string filename, Strain& strain)
{
    set_regional_area_strain(strain); 
    std::ofstream datafile;
    datafile.open(std::filesystem::path(DATA_DIR) / filename, std::ios_base::app);
    
    // Write the header 
    datafile << "\nGlobal Area Strain, Inferior Free Wall Area Strain, Lateral Free Wall Area Strain, "
             << "Anterior Free Wall Area Strain, Septal Body Area Strain \n"; 
   
    // Find global area strain 
    strain.global_area_strain = area_strain(strain.global_ed_area, strain.global_es_area); 
    datafile << strain.global_area_strain << ", "; 

    // Display area strain for each region  
    for (auto& a : strain.strain_values) {
        datafile << a << ", ";
    }

     datafile << "\n\nGlobal Longitudinal Strain, Inferior Free Wall Longitudinal Strain, Lateral Free Wall Longitudinal Strain, Anterior Free Wall Longitudinal Strain, Septal Body Longitudinal Strain\n"; 

    // Find global and regional longitudinal strain 
    std::string meshFile = construct_file_string(1); 
    Mesh rv = get_mesh(meshFile); 
    strain.longitudinal_strain = longitudinal_strain(strain.vertices_es, strain.vertices_ed, rv.triangles, rv.vertexToTri, strain.long_axis, strain); 

    datafile << strain.longitudinal_strain << ", "; // global longitudinal strain
    for (auto& la : strain.l_strain_values) {
        datafile << la << ", "; 
    }
    
    datafile << "\n\nGlobal Radial Strain, Inferior Free Wall Radial Strain, Lateral Free Wall Radial Strain, Anterior Free Wall Radial Strain, Septal Body Radial Strain\n"; 
    strain.radial_axis = 1.0f * ((float)rv.rad_radius) * glm::normalize(find_radial_axis(strain.vertices_ed, strain.long_axis)); 
    strain.radial_strain = radial_strain(strain.vertices_es, strain.vertices_ed, strain); 

    datafile << strain.radial_strain << ", "; // global radial strain
    for (auto& la : strain.r_strain_values) {
        datafile << la << ", ";
    }

    datafile << "\n\nGlobal Circumferential Strain, Inferior Free Wall Circumferential Strain, Lateral Free Wall Circumferential Strain, Anterior Free Wall Circumferential Strain, Septal Body Circumferential Strain\n"; 
    strain.circ_axis = 1.0f * ((float)rv.circ_radius) * glm::normalize(find_circumferential_axis(strain.vertices_ed, strain.long_axis, strain.radial_axis));  
    strain.circumferential_strain = circumferential_strain(strain.vertices_es, strain.vertices_ed, strain); 

    datafile << strain.circumferential_strain << ", "; // global circumferential strain
    for (auto& la : strain.c_strain_values) {
        datafile << la << ", ";
    }
          
    std::cout << "radial axis: " << strain.radial_axis.x << " " << strain.radial_axis.y << " " << strain.radial_axis.z << std::endl;
    std::cout << "circ axis: " << strain.circ_axis.x << " " << strain.circ_axis.y << " " << strain.circ_axis.z << std::endl;
    std::cout << "dot product c and r: " << glm::dot(strain.circ_axis, strain.radial_axis) << std::endl;
    std::cout << "dot product r and l: " << glm::dot(strain.radial_axis, strain.long_axis) << std::endl;
    std::cout << "dot product c and l: " << glm::dot(strain.circ_axis, strain.long_axis) << std::endl; 

    datafile.close();
}

/*
* This main function for calculations
*/
int main_calculations(std::string csvfile, int frames)
{
    // File name strings 
    std::string fileName;
    
    // Variables used for calculating strain 
    double min_vol = 10000; // end systole
    double max_vol = 0;  // end diastole
    Strain strain = {}; 

    // Apex point 
    int lower_point = 906;
    // Center of inflow tract
    int c1 = 102; 
    // Center of outflow tract 
    int c2 = 63;

    for (int i = 0; i <= frames; i++) {      
        fileName = construct_file_string(i); 
        Mesh rv = get_mesh(fileName); 

        // Printing Info for CSV file
        PrintInfo print_info = {}; 

        // Calculate global quantities 
        print_info.volume = find_volume(rv.triangles, rv.vertices) / 1000;
        print_info.surface_area = find_surface_area(rv.triangles, rv.vertices);
        print_info.curvature = find_curvature(rv.triangles, rv.vertices, rv.vertexToTri);
        print_info.index_curv = find_indexed_curvature(rv.vertices); 

        // Calculate regional quantities 
        print_info.curvatures = find_regional_curvature(rv.vertices);
        print_info.volumes = regional_volumes(rv.vertices, rv.triangles, rv.vertexToTri); 
        print_info.surface_areas = regional_surface_areas(rv.vertices, rv.triangles, rv.vertexToTri); 

        // Write to CSV file
        print_info.minmax = find_min_max(rv.vertices); 
        write_to_file(csvfile, print_info, i);

        // Recording variables for calculating strain
        if (print_info.volume < min_vol) { // end-systole, min. volume
            strain.global_es_area = print_info.surface_area; 
            min_vol = print_info.volume; 
            strain.es_areas = std::vector<double>(print_info.surface_areas.begin(), print_info.surface_areas.end()); 
            strain.vertices_es = std::vector<Vertex>(rv.vertices.begin(), rv.vertices.end()); 
        } else if (print_info.volume > max_vol) { // end diastole, max. volume
            strain.global_ed_area = print_info.surface_area; 
            max_vol = print_info.volume; 
            strain.ed_areas = std::vector<double>(print_info.surface_areas.begin(), print_info.surface_areas.end()); 
            strain.vertices_ed = std::vector<Vertex>(rv.vertices.begin(), rv.vertices.end());       
            strain.long_axis = find_long_axis(rv.vertices, lower_point, c1, c2); 
        }
    }
    std::cout << "here" << std::endl;
    // Output and _sets_ strain calculations 
    write_strain_to_file(csvfile, strain); 

    return 0; 
}

/*
 * This main function would be to plot the RV beutel
 */
int main_visual()
{
    Window window { "RV Beutel Visualisation", glm::ivec2(1000), OpenGLVersion::GL2 };
    std::string fileName = "healthy-1/Young healthy volunteer_000.obj";

    std::string ring = "ring-indices.txt"; // ring-indices, ring-sphere ring-large
    std::string exclude_vertices = "exclude.txt";
    std::string regions = "region-v2.txt";

    Trackball trackball { &window, glm::radians(60.0f), 2.0f, 0.387463093f, -0.293215364f };
    trackball.disableTranslation();
    printHelp();

    // Load the mesh file and ring file
    std::ifstream ifile;
    ifile.open(std::filesystem::path(DATA_DIR) / fileName);
    Mesh rv = loadMeshRV(ifile);
    loadRingFromFile(ring, rv.vertices);

    // ImGUI lights
    ProgramState state {};
    state.myMesh = rv;
    // Extra processing needed for RV Beutel
    if (rv.vertices.size() > 937) {
        // Vertices used for calculations 
        center_mesh(rv.vertices);
        mark_regions(regions, rv.vertices);
        mark_excluded(exclude_vertices, rv.vertices);
        // Vertices used for GUI 
        scale_mesh(state.myMesh.vertices);
        center_mesh(state.myMesh.vertices);
        mark_regions(regions, state.myMesh.vertices);
        mark_excluded(exclude_vertices, state.myMesh.vertices); 
    }

    // Create background and lights for mesh 
    state.materialInformation.Kd = glm::vec3(75, 139, 59) / 255.0f;
    state.materialInformation.Ks = glm::vec3(221, 42, 116) / 255.0f;
    state.materialInformation.shininess = 20.0f;
    meshFlipZ(state.myMesh);
    float lp = .85f;
    state.lights.push_back(Light { glm::vec3(lp, lp, lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(-lp, lp, lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(lp, -lp, lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(lp, lp, -lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(-lp, -lp, lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(lp, -lp, -lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(-lp, lp, -lp), glm::vec3(224, 215, 73) / 255.0f });
    state.lights.push_back(Light { glm::vec3(-lp, -lp, -lp), glm::vec3(224, 215, 73) / 255.0f });

    // Calculate the actual volume captured by mesh
    RVInfo info {};
    info.volume = find_volume(rv.triangles, rv.vertices) / 1000;
    info.surfaceArea = find_surface_area(rv.triangles, rv.vertices);
    info.curvature = find_curvature(rv.triangles, rv.vertices, rv.vertexToTri);

    // Reset to show _indexed_ curvature 
    info.curvature = find_indexed_curvature(rv.vertices); 
    info.radius = rv.circ_radius;

    // GUI (Visual Debug) features
    set_regional(info, rv.vertices);
    std::vector<glm::vec3> vertexColors = heat_color(rv.triangles, rv.vertices, rv.vertexToTri);
    std::vector<Ray> normals = find_normals(state.myMesh.vertices, state.myMesh.triangles);

    // Window Pop-up
    window.registerCharCallback([&](unsigned unicodeCodePoint) {
        keyboard(static_cast<unsigned char>(unicodeCodePoint), state, window, trackball);
    });

    while (!window.shouldClose()) {
        window.updateInput();
        glViewport(0, 0, window.getFrameBufferSize().x, window.getFrameBufferSize().y);

        std::vector<glm::vec3> lightColors(state.myMesh.vertices.size());
        computeLighting(state, trackball.position(), lightColors);

        draw(state, trackball, vertexColors, state.myMesh.vertices);
        drawUI(state, trackball, info, normals);

        window.swapBuffers();
    }

    return 0; 
}

// Main function 
int main(int argc, char** argv)
{
    // Either main_calculations or main_visual
    //main_visual(); 
    main_calculations("data-3.csv", 38);
}