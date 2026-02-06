# Add the RunningGate library
cmake_minimum_required(VERSION 3.10)
project(RunningGate VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-g -fmax-errors=5 -O3)


# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core Hist Graf Gpad)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Add subdirectories
add_subdirectory(modules/theuerkauf_fitter)
add_subdirectory(src)

# Add the SampleCode executable
add_executable(running_gate ExampleCode.cpp)

# Include the src directory so SampleCode can find RunningGate.h
target_include_directories(running_gate PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Link the executable with RunningGate library and ROOT
target_link_libraries(running_gate PRIVATE RunningGate ${ROOT_LIBRARIES})

# Create symbolic link to example_data directory in build directory
add_custom_target(create_example_data_symlink ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_SOURCE_DIR}/example_data
        ${CMAKE_BINARY_DIR}/example_data
    COMMENT "Creating symbolic link to example_data directory"
)
