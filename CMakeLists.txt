cmake_minimum_required(VERSION 3.18)
set (CMAKE_CXX_STANDARD 17)

project(bcd)

if(${WITH_QT})
    find_package(CGAL COMPONENTS Core Qt6 REQUIRED) # test: CGAL 5.6.1
else()
    find_package(CGAL COMPONENTS Core REQUIRED) # test: CGAL 5.6.1
endif()
find_package(Eigen3 REQUIRED)

find_package(RapidJSON)
find_package(PCL 1.14 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
find_package(NLopt HINTS ${CMAKE_INSTALL_PREFIX}/lib/cmake)
find_package(sdf)

if (APPLE)
    link_directories(/usr/local/lib)
    find_package(LLVM REQUIRED CONFIG)
    set(CMAKE_C_COMPILER "${LLVM_TOOLS_BINARY_DIR}/clang")
    set(CMAKE_CXX_COMPILER "${LLVM_TOOLS_BINARY_DIR}/clang++")
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    find_package(OpenMP REQUIRED)
elseif(UNIX)
    find_package(OpenMP)
endif()

add_executable(bcd 
    main.cpp util.hpp global_variables.hpp reconstruction.hpp specific_search.hpp 
    bnb.hpp Node.hpp Component.hpp Edge.hpp)

if(${DEBUG_SEARCH})
    target_compile_definitions(bcd PRIVATE DEBUG_SEARCH=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_SEARCH=0)
endif()

if(${DEBUG_CUT})
    target_compile_definitions(bcd PRIVATE DEBUG_CUT=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_CUT=0)
endif()

if(${DEBUG_NODE})
    target_compile_definitions(bcd PRIVATE DEBUG_NODE=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_NODE=0)
endif()

if(${DEBUG_CONCAVE_EDGE})
    target_compile_definitions(bcd PRIVATE DEBUG_CONCAVE_EDGE=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_CONCAVE_EDGE=0)
endif()

if(${DEBUG_SYMMETRY})
    target_compile_definitions(bcd PRIVATE DEBUG_SYMMETRY=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_SYMMETRY=0)
endif()

if(${DEBUG_CENTRAL})
    target_compile_definitions(bcd PRIVATE DEBUG_CENTRAL=1)
else()
    target_compile_definitions(bcd PRIVATE DEBUG_CENTRAL=0)
endif()

if(${LOG_TIME})
    target_compile_definitions(bcd PRIVATE LOG_TIME=1)
else()
    target_compile_definitions(bcd PRIVATE LOG_TIME=0)
endif()


if(${LOG_SEARCH})
    target_compile_definitions(bcd PRIVATE LOG_SEARCH=1)
else()
    target_compile_definitions(bcd PRIVATE LOG_SEARCH=0)
endif()

if(${VIS_CONCAVE_CUT})
    target_compile_definitions(bcd PRIVATE VIS_CONCAVE_CUT=1)
else()
    target_compile_definitions(bcd PRIVATE VIS_CONCAVE_CUT=0)
endif()

if(${WITH_QT})
    target_compile_definitions(bcd PRIVATE VIS=1)    

    if(${VIS_NODE})
        target_compile_definitions(bcd PRIVATE VIS_NODE=1)
    else()
        target_compile_definitions(bcd PRIVATE VIS_NODE=0)
    endif()

    if(${VIS_CHILD})
        target_compile_definitions(bcd PRIVATE VIS_CHILD=1)
    else()
        target_compile_definitions(bcd PRIVATE VIS_CHILD=0)
    endif()

    if(${VIS_LOD})
        target_compile_definitions(bcd PRIVATE VIS_LOD=1)
    else()
        target_compile_definitions(bcd PRIVATE VIS_LOD=0)
    endif()

    target_link_libraries(bcd
                CGAL::CGAL CGAL::CGAL_Qt6
                OpenMP::OpenMP_CXX 
                sdf::sdf
                ${PCL_LIBRARIES}
                NLopt::nlopt)
else()
    target_link_libraries(bcd
                CGAL::CGAL OpenMP::OpenMP_CXX 
                sdf::sdf
                ${PCL_LIBRARIES}
                NLopt::nlopt)
endif()