cmake_minimum_required(VERSION 3.14)

# Name the project
project(ProgrammingProject)

# Specify the c++ standard
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 14)

# Point to the /src/ subdirectory
#add_subdirectory(src)

# Find needed Lapack package
find_package(LAPACK REQUIRED)

# Add the main executable
add_executable(Solution ./src/molecule.cpp ./src/Solution.cpp)

# Link it to LAPACK for matrix operations
target_link_libraries(Solution ${LAPACK_LIBRARIES})

# Add HF exeutable
add_executable(HF ./src/hf.cpp ./src/hf_wfn.cpp ./src/tei_handler.cpp)
