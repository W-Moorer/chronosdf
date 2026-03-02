cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0091 NEW)

project(chrono_sdf_contact LANGUAGES CXX)

option(SDF_CONTACT_BUILD_DEMOS "Build SDF contact skeleton demos" ON)

find_package(Chrono CONFIG REQUIRED OPTIONAL_COMPONENTS FEA)

if(NOT Chrono_FOUND)
    message(FATAL_ERROR "Could not find Chrono")
endif()

add_library(chrono_sdf_contact
    src/SdfVolume.cpp
    src/SdfOctreeVolume.cpp
    src/SdfRegistry.cpp
    src/SdfPairCache.cpp
    src/SdfNarrowphase.cpp
    src/SdfClustering.cpp
    src/SdfNarrowphaseCallback.cpp
    src/SdfCustomCollisionCallback.cpp
)

target_compile_features(chrono_sdf_contact PUBLIC cxx_std_17)
target_include_directories(chrono_sdf_contact PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(chrono_sdf_contact PUBLIC ${CHRONO_TARGETS})

if(MSVC)
    set_target_properties(chrono_sdf_contact PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    target_compile_options(chrono_sdf_contact PRIVATE /Zc:inline)
endif()

if(SDF_CONTACT_BUILD_DEMOS)
    add_executable(demo_MBS_sdf_contact demos/demo_MBS_sdf_contact.cpp)
    target_link_libraries(demo_MBS_sdf_contact PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_contact PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_contact PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    if(CHRONO_FEA_FOUND)
        add_executable(demo_FEA_sdf_mesh_contact demos/demo_FEA_sdf_mesh_contact.cpp)
        target_link_libraries(demo_FEA_sdf_mesh_contact PRIVATE chrono_sdf_contact)
        target_compile_definitions(demo_FEA_sdf_mesh_contact PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
        if(MSVC)
            set_target_properties(demo_FEA_sdf_mesh_contact PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
        endif()
    else()
        message(STATUS "Chrono::FEA not found; skipping demo_FEA_sdf_mesh_contact")
    endif()

    add_executable(demo_MBS_sdf_trimesh_contact demos/demo_MBS_sdf_trimesh_contact.cpp)
    target_link_libraries(demo_MBS_sdf_trimesh_contact PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_trimesh_contact PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_trimesh_contact PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_sdf_contact demos/demo_MBS_sdf_sdf_contact.cpp)
    target_link_libraries(demo_MBS_sdf_sdf_contact PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_sdf_contact PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_sdf_contact PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_mesh_baseline_compare demos/demo_MBS_sdf_mesh_baseline_compare.cpp)
    target_link_libraries(demo_MBS_sdf_mesh_baseline_compare PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_mesh_baseline_compare PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_mesh_baseline_compare PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_sdf_stability_compare demos/demo_MBS_sdf_sdf_stability_compare.cpp)
    target_link_libraries(demo_MBS_sdf_sdf_stability_compare PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_sdf_stability_compare PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_sdf_stability_compare PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_octree_pareto demos/demo_MBS_sdf_octree_pareto.cpp)
    target_link_libraries(demo_MBS_sdf_octree_pareto PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_octree_pareto PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_octree_pareto PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_p2_scene_matrix demos/demo_MBS_sdf_p2_scene_matrix.cpp)
    target_link_libraries(demo_MBS_sdf_p2_scene_matrix PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_p2_scene_matrix PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_p2_scene_matrix PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    add_executable(demo_MBS_sdf_p2_scaleup demos/demo_MBS_sdf_p2_scaleup.cpp)
    target_link_libraries(demo_MBS_sdf_p2_scaleup PRIVATE chrono_sdf_contact)
    target_compile_definitions(demo_MBS_sdf_p2_scaleup PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    if(MSVC)
        set_target_properties(demo_MBS_sdf_p2_scaleup PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()

    if(WIN32)
        add_CHRONO_DLLS_copy_command()
        if(TARGET COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_contact COPY_CHRONO_DLLS)
            if(TARGET demo_FEA_sdf_mesh_contact)
                add_dependencies(demo_FEA_sdf_mesh_contact COPY_CHRONO_DLLS)
            endif()
            add_dependencies(demo_MBS_sdf_trimesh_contact COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_sdf_contact COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_mesh_baseline_compare COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_sdf_stability_compare COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_octree_pareto COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_p2_scene_matrix COPY_CHRONO_DLLS)
            add_dependencies(demo_MBS_sdf_p2_scaleup COPY_CHRONO_DLLS)
        endif()
    endif()
endif()
