cmake_minimum_required(VERSION 3.0)

# name of application. replace 'app' with desired app name
set(APP_NAME casm_viewer)

project(${APP_NAME})

add_subdirectory(external/tinc)

# CASM viewer application
add_executable(${APP_NAME}
  src/main.cpp
  src/datadisplay.cpp
  src/datasetmanager.cpp
  src/modalsynth.cpp

  src/datadisplay.hpp
  src/datasetmanager.hpp
  src/modalsynth.hpp
  )

target_link_libraries(${APP_NAME} PRIVATE tinc)

#if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/al_ext)
#  message("Buiding extensions in al_ext")
#  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/al_ext)
#  get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)
#  target_link_libraries(${APP_NAME} PRIVATE ${AL_EXT_LIBRARIES})
#endif()

# compile flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
set(app_compile_flags
    /bigobj /F2000000
)
endif()

# binaries are put into the ./bin directory by default
set_target_properties(${APP_NAME} PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/bin
)

# --- Tools ########################################
add_executable(ECI_explorer
  src/eci_explorer.cpp
  )

target_link_libraries(ECI_explorer PRIVATE tinc)

# binaries are put into the ./bin directory by default
set_target_properties(ECI_explorer PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/bin
)

