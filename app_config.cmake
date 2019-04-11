# allolib location
set(allolib_directory allolib)

# list your app files here
set(app_files_list
  vdv_group/al_DataScript.hpp
  vdv_group/common.hpp
  vdv_group/processors.hpp
  vdv_group/parameterspace.hpp
  vdv_group/datasetmanager.hpp
  vdv_group/instanced_mesh.hpp
  vdv_group/simulator.cpp
  vdv_group/al_VASPReader.hpp
  vdv_group/json.hpp
  vdv_group/slice.hpp
)

# other directories to include
set(app_include_dirs
  ${CMAKE_CURRENT_LIST_DIR}
  ${allolib_directory}/external/Gamma
)

# other libraries to link
set(app_link_libs
)

# definitions
set(app_definitions
)

# compile flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
set(app_compile_flags
    /bigobj /F2000000
)
endif()

# linker flags, with `-` in the beginning
set(app_linker_flags
)
