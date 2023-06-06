cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(ComputerGraphics C CXX)

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/framework")
	# Create framework library and include CMake scripts (compiler warnings, sanitizers and static analyzers).
	add_subdirectory("framework")
else()
	# During development the framework lives in parent folder.
	add_subdirectory("../../../framework/" "${CMAKE_BINARY_DIR}/framework/")
endif()
find_package(OpenGL REQUIRED)

add_executable(BeutelAnalysis
	"src/main.cpp"
	"src/draw.cpp" "src/mesh_loader.h" "src/mesh_loader.cpp" "src/util.h" "src/util.cpp" 
	"src/strain_calculation.h" "src/strain_calculation.cpp")
target_compile_features(BeutelAnalysis PRIVATE cxx_std_20)
target_link_libraries(BeutelAnalysis PRIVATE CGFramework OpenGL::GLU)
enable_sanitizers(BeutelAnalysis)
set_project_warnings(BeutelAnalysis)

target_compile_definitions(BeutelAnalysis PRIVATE "-DDATA_DIR=\"${CMAKE_CURRENT_LIST_DIR}/data/\"")

