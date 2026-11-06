add_library(ccell_preprocess)
add_library(ccell::preprocess ALIAS ccell_preprocess)
set_target_properties(ccell_preprocess PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  EXPORT_NAME preprocess_binary)
target_compile_features(ccell_preprocess PUBLIC cxx_std_20)
target_sources(ccell_preprocess
  PUBLIC FILE_SET CXX_MODULES BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src"
    FILES src/preprocess.cppm
  PRIVATE src/preprocess.cpp src/hvg.cpp)
target_link_libraries(ccell_preprocess PUBLIC anndatapp::anndatapp)
target_compile_options(ccell_preprocess PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS ccell_preprocess EXPORT ccellTargets
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_DATADIR}/ccell/preprocess/modules
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(CCELL_BUILD_TESTS)
  foreach(suite IN ITEMS numerics hvg)
    add_executable(ccell_preprocess_${suite} tests/test_${suite}.cpp)
    target_link_libraries(ccell_preprocess_${suite} PRIVATE ccell::preprocess)
    add_test(NAME ccell-preprocess-${suite} COMMAND ccell_preprocess_${suite})
  endforeach()

  add_executable(ccell_preprocess_example examples/basic.cpp)
  target_link_libraries(ccell_preprocess_example PRIVATE ccell::preprocess)
  add_test(NAME ccell-preprocess-example COMMAND ccell_preprocess_example)
endif()

if(CCELL_BUILD_QUALIFICATION)
  add_executable(ccell_preprocess_interop_driver tests/interop_driver.cpp)
  target_link_libraries(ccell_preprocess_interop_driver PRIVATE ccell::preprocess)
  add_executable(ccell_preprocess_benchmark benchmarks/qualification.cpp)
  target_link_libraries(ccell_preprocess_benchmark PRIVATE ccell::preprocess)
endif()
