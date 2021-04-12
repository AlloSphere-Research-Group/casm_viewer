cmake_minimum_required(VERSION 3.0)

# name of application. replace 'app' with desired app name
set(APP_NAME casm_viewer)

project(${APP_NAME})

# path to main source file
add_executable(${APP_NAME}
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

# add allolib as a subdirectory to the project
add_subdirectory(external/tinc)

# link allolib to project
#target_link_libraries(${APP_NAME} PRIVATE al)
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

# example line for find_package usage
# find_package(Qt5Core REQUIRED CONFIG PATHS "C:/Qt/5.12.0/msvc2017_64/lib" NO_DEFAULT_PATH)

# replace ${PATH_TO_INCLUDE_DIR} before including other libraries
# target_include_directories(${APP_NAME} PRIVATE ${PATH_TO_INCLUDE_DIR})

# replace ${PATH_TO_LIB_FILE} before linking other libraries
# target_link_libraries(${APP_NAME} PRIVATE ${PATH_TO_LIB_FILE})

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


option(BUILD_CASMVIEWER_RELEASE "" OFF)
if (BUILD_CASMVIEWER_RELEASE)

  if(WIN32)
    set(WIN_RES res/myrioi-windows.rc)
  endif(WIN32)

  # ----- Packaging
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX  ${CMAKE_CURRENT_LIST_DIR}/package CACHE PATH "..." FORCE)
  endif()


  if(WIN32 AND NOT UNIX)
#      install(TARGETS myrioi-vr
#        BUNDLE DESTINATION . COMPONENT Runtime
#        RUNTIME DESTINATION bin COMPONENT Runtime
#        )
#    install(CODE [[
#        file(GET_RUNTIME_DEPENDENCIES
#            RESOLVED_DEPENDENCIES_VAR RES
#            UNRESOLVED_DEPENDENCIES_VAR UNRES
#            CONFLICTING_DEPENDENCIES_PREFIX CONFLICTING_DEPENDENCIES
#            EXECUTABLES myrioi-vr
#        )]]
#    )

#    if (EXISTS "C:/Program Files/Mega-Nerd/libsndfile/bin/libsndfile-1.dll")
#        add_custom_command(
#                TARGET myrioi-vr POST_BUILD
#                COMMAND ${CMAKE_COMMAND} -E copy
#                        "C:/Program Files/Mega-Nerd/libsndfile/bin/libsndfile-1.dll"
#                        "${CMAKE_CURRENT_LIST_DIR}/myrioi-vr")
#    endif()

#    get_target_property(OPENVR_ROOT al_openvr OPENVR_ROOT)
#    if (EXISTS "${OPENVR_ROOT}/bin/win64/openvr_api.dll")
#        add_custom_command(
#                TARGET myrioi-vr POST_BUILD
#                COMMAND ${CMAKE_COMMAND} -E copy
#                        "${OPENVR_ROOT}/bin/win64/openvr_api.dll"
#                        "${CMAKE_CURRENT_LIST_DIR}/myrioi-vr")
#    endif()

  elseif(APPLE)

    # Note Mac specific extension .app
    set(APPS "${CMAKE_CURRENT_LIST_DIR}/CASM_Viewer.app")

    # Directories to look for dependencies
    set(DIRS ${CMAKE_BINARY_DIR})

    install(CODE "include(BundleUtilities)
      fixup_bundle(\"${APPS}\" \"\" \"${DIRS}\")")

    set(CPACK_GENERATOR "DragNDrop")

    include(InstallRequiredSystemLibraries)

    set(CPACK_BUNDLE_NAME "CASM_Viewer")
    set(CPACK_BUNDLE_PLIST "${CMAKE_CURRENT_LIST_DIR}/res/macos/Info.plist")
    set(CPACK_BUNDLE_ICON "${CMAKE_CURRENT_LIST_DIR}/res/macos/myrioi_icon.icns")
    set(CPACK_PACKAGE_FILE_NAME "CASM_Viewer")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CASM_Viewer")
    set(CPACK_PACKAGE_VENDOR "Allosphere Research Group, Van der Ven Research group")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_LIST_DIR}/res/readme.md")
    set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_LIST_DIR}/res/macos/myrioi_icon.icns")
    #  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
    set(CPACK_PACKAGE_VERSION_MAJOR "1")
    set(CPACK_PACKAGE_VERSION_MINOR "0")
    set(CPACK_PACKAGE_VERSION_PATCH "0")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

    set(CPACK_PACKAGE_EXECUTABLES "casm_viewer")

    include(CPack)

  else() # Linux

  endif()


endif(BUILD_CASMVIEWER_RELEASE)

