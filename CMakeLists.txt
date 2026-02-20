cmake_minimum_required(VERSION 3.10)
project(OpenVort C)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find libconfig using PkgConfig
find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBCONFIG REQUIRED libconfig)

# Compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -Wall -Wextra -pedantic -ggdb")

# Include directories
include_directories(include ${LIBCONFIG_INCLUDE_DIRS})

# Source files
set(SOURCES
    src/boundary_images.c
    src/configuration.c
    src/external_velocity.c
    src/initial_conditions.c
    src/main.c
    src/octree.c
    src/tangle.c
    src/util.c
    src/vec3_maths.c
    src/vortex_constants.c
    src/vortex_dynamics.c
    src/vortex_injection.c
    src/vortex_utils.c
)

# Executable
add_executable(vortices ${SOURCES})

# Compile definitions
target_compile_definitions(vortices PRIVATE _DEBUG_)

# Link libraries
target_link_libraries(vortices OpenMP::OpenMP_C m ${LIBCONFIG_LIBRARIES})
