cmake_minimum_required(VERSION 3.28)
project(ccell VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(CCELL_BUILD_TESTS "Build module composition and installed-consumer checks" ${PROJECT_IS_TOP_LEVEL})
option(CCELL_BUILD_QUALIFICATION "Build preprocessing numerical qualification drivers" OFF)
if(CCELL_BUILD_TESTS)
  enable_testing()
endif()

add_subdirectory(anndata)
add_subdirectory(preprocess)
add_subdirectory(velocity)

add_library(ccell)
add_library(ccell::ccell ALIAS ccell)
set_target_properties(ccell PROPERTIES POSITION_INDEPENDENT_CODE ON EXPORT_NAME ccell_binary)
target_compile_features(ccell PUBLIC cxx_std_20)
target_sources(ccell PUBLIC FILE_SET CXX_MODULES BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}"
  FILES ccell.ccm src/ccell.anndata.cppm src/ccell.preprocess.cppm src/ccell.velocity.cppm)
foreach(dependency IN ITEMS anndatapp::anndatapp ccell::preprocess velocitypp::velocitypp)
  target_link_libraries(ccell PUBLIC
    $<BUILD_INTERFACE:${dependency}> $<INSTALL_INTERFACE:${dependency}>)
endforeach()

install(TARGETS ccell EXPORT ccellTargets
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_DATADIR}/ccell/modules
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT ccellTargets NAMESPACE ccell:: CXX_MODULES_DIRECTORY cxx-modules
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ccell)
configure_package_config_file(cmake/ccellConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ccellConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ccell)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/ccellConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ccellConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ccellConfigVersion.cmake" cmake/ccellModuleProvider.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ccell)
install(FILES LICENSE THIRD_PARTY_NOTICES.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/ccell)

if(CCELL_BUILD_TESTS)
  foreach(suite IN ITEMS umbrella components)
    add_executable(ccell_${suite} tests/${suite}.cpp)
    target_link_libraries(ccell_${suite} PRIVATE ccell::ccell)
    add_test(NAME ccell-${suite} COMMAND ccell_${suite})
  endforeach()
  add_test(NAME ccell-installed-consumer COMMAND ${CMAKE_COMMAND}
    -DCCELL_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DCCELL_BUILD_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DCCELL_COMPILER=${CMAKE_CXX_COMPILER}
    -DCCELL_CONFIG=$<CONFIG>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test-installed-consumer.cmake)
endif()
