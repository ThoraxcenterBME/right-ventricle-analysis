#include "draw.h"
#include "util.h"
#include "shade.h"
#include "mesh_loader.h"
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

struct RVInfo {
    float volume;
    float surfaceArea; 
    float curvature; 
    float radius; 
    float region1; 
    float region2; 
    float region3; 
    float region4; 
    float region5; 
    float region6; 
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

void draw(const ProgramState& state, const Trackball& camera, std::vector<glm::vec3> vertexColors)
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
    drawAxis();

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
    std::string curveString = "Curvature: " + std::to_string(info.curvature);
    ImGui::Text(curveString.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Check for drawing normals
    ImGui::Checkbox("Show Surface Normals", &show_normal_gui); 
    if (show_normal_gui) {
        draw_normal_rays(normals); 
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
    std::string r1_curvature = "Curvature (Inflow Tract): " + std::to_string(info.region1);
    ImGui::TextColored(ImVec4(color.c1.x, color.c1.y, color.c1.z, 1.0), r1_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 2
    std::string r2_curvature = "Curvature (Outflow Tract): " + std::to_string(info.region2);
    ImGui::TextColored(ImVec4(color.c2.x, color.c2.y, color.c2.z, 1.0), r2_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 3
    std::string r3_curvature = "Curvature (Septal Body): " + std::to_string(info.region3);
    ImGui::TextColored(ImVec4(color.c3.x, color.c3.y, color.c3.z, 1.0), r3_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();


    // Display Region 4
    std::string r4_curvature = "Curvature (Free-wall Body): " + std::to_string(info.region4);
    ImGui::TextColored(ImVec4(color.c4.x, color.c4.y, color.c4.z, 1.0), r4_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 5
    std::string r5_curvature = "Curvature (Septal Apex): " + std::to_string(info.region5);
    ImGui::TextColored(ImVec4(color.c5.x, color.c5.y, color.c5.z, 1.0), r5_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display Region 6
    std::string r6_curvature = "Curvature (Free-wall Apex): " + std::to_string(info.region6);
    ImGui::TextColored(ImVec4(color.c6.x, color.c6.y, color.c6.z, 1.0), r6_curvature.c_str());
    ImGui::Spacing();
    ImGui::Separator();

    // Display other information ...
    ImGui::End();
}

void printHelp()
{
    std::cout << "Program Usage:" << std::endl;
    std::cout << "SPACE - replaces mouse click for selection, will then call your light placement function" << std::endl;
}

void set_regional(RVInfo& info, std::vector<Vertex>& vertices) {
    std::vector<double> curvatures = find_regional_curvature(vertices); 
    info.region1 = curvatures[0];
    info.region2 = curvatures[1];
    info.region3 = curvatures[2];
    info.region4 = curvatures[3];
    info.region5 = curvatures[4];
    info.region6 = curvatures[5];
}

/*
* This main function would be to plot the RV beutel 
*/
int main(int argc, char** argv)
{
    Window window { "RV Beutel Visualisation", glm::ivec2(1000), OpenGLVersion::GL2 };
    std::string fileName = "ref.obj";
    std::string ring = "ring-indices.txt"; // ring-indices, ring-sphere ring-large
    std::string exclude_vertices = "exclude.txt"; 
    std::string regions = "region.txt"; 

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
        center_mesh(rv.vertices);
        mark_excluded(exclude_vertices, rv.vertices);
        mark_regions(regions, rv.vertices);
        scale_mesh(state.myMesh.vertices);
        center_mesh(state.myMesh.vertices);
        mark_regions(regions, state.myMesh.vertices);
    }

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
    info.volume = find_volume(rv.triangles, rv.vertices);
    info.surfaceArea = find_surface_area(rv.triangles, rv.vertices); 
    info.curvature = find_curvature(rv.triangles, rv.vertices, rv.vertexToTri);
    info.radius = rv.radius; 

    // GUI (Visual Debug) features 
    std::vector<glm::vec3> vertexColors = heat_color(rv.triangles, rv.vertices, rv.vertexToTri);
    set_regional(info, rv.vertices);
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

        draw(state, trackball, vertexColors);
        drawUI(state, trackball, info, normals);
       
        window.swapBuffers();
    }
}
