cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(ComputerGraphics C CXX)

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/framework")
	# Create framework library and include CMake scripts (compiler warnings, sanitizers and static analyzers).
	add_subdirectory("framework")
else()
	# During development the framework lives in parent folder.
	add_subdirectory("../../../framework/" "${CMAKE_BINARY_DIR}/framework/")
endif()

add_executable(ShadingAssignment
	"src/main.cpp"
	"src/draw.cpp" "src/mesh_loader.h" "src/mesh_loader.cpp" "src/util.h" "src/util.cpp")
target_compile_features(ShadingAssignment PRIVATE cxx_std_20)
target_link_libraries(ShadingAssignment PRIVATE CGFramework)
enable_sanitizers(ShadingAssignment)
set_project_warnings(ShadingAssignment)

target_compile_definitions(ShadingAssignment PRIVATE "-DDATA_DIR=\"${CMAKE_CURRENT_LIST_DIR}/data/\"")

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/grading_tests/")
	add_subdirectory("grading_tests")
endif()

