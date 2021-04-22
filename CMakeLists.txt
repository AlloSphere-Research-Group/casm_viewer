cmake_minimum_required(VERSION 3.0)

# name of application. replace 'app' with desired app name
set(APP_NAME casm_viewer)

project(${APP_NAME})

# path to main source file
add_executable(${APP_NAME} MACOSX_BUNDLE
  src/main.cpp
  src/datadisplay.cpp
  src/datasetmanager.cpp
  src/slice.cpp
  src/modalsynth.cpp

  src/datadisplay.hpp
  src/datasetmanager.hpp
  src/slice.hpp
  src/modalsynth.hpp
  )

# add tinc as a subdirectory to the project
add_subdirectory(external/tinc)

# link tinc to project
target_link_libraries(${APP_NAME} PRIVATE tinc)

# binaries are put into the ./bin directory by default
set_target_properties(${APP_NAME} PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/bin
)

################################################
# Deployment and packaging

# ----- Packaging
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX  ${CMAKE_CURRENT_LIST_DIR}/package CACHE PATH "..." FORCE)
endif()

install(TARGETS casm_viewer
  BUNDLE DESTINATION . COMPONENT Runtime
  RUNTIME DESTINATION bin COMPONENT Runtime
  )

# Note Mac specific extension .app
set(APPS "${CMAKE_CURRENT_LIST_DIR}/bin/casm_viewer.app")

# Directories to look for dependencies
set(DIRS ${CMAKE_BINARY_DIR})

install(CODE "include(BundleUtilities)
        set(BU_CHMOD_BUNDLE_ITEMS TRUE)
        fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")")

set(CPACK_GENERATOR "DragNDrop")

include(InstallRequiredSystemLibraries)

set(CPACK_BUNDLE_NAME "CASM_Viewer")
set(CPACK_BUNDLE_PLIST "${CMAKE_CURRENT_LIST_DIR}/res/macos/Info.plist")
set(CPACK_BUNDLE_ICON "${CMAKE_CURRENT_LIST_DIR}/res/macos/myrioi_icon.icns")
set(CPACK_PACKAGE_FILE_NAME "CASM_Viewer")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CASM_Viewer")
set(CPACK_PACKAGE_VENDOR "Allosphere Research Group, Van der Ven Research group")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_LIST_DIR}/readme.md")
set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_LIST_DIR}/res/macos/myrioi_icon.icns")
#  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

set(CPACK_PACKAGE_EXECUTABLES "casm_viewer")

include(CPack)