add_library(ccell_preprocess)
add_library(ccell::preprocess ALIAS ccell_preprocess)
set_target_properties(ccell_preprocess PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  EXPORT_NAME preprocess_binary)
target_compile_features(ccell_preprocess PUBLIC cxx_std_20)
target_sources(ccell_preprocess
  PUBLIC FILE_SET CXX_MODULES BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}"
    FILES preprocess.cppm qc.cppm hvg.cppm
  PRIVATE preprocess.cpp qc.cpp hvg.cpp)
target_link_libraries(ccell_preprocess PUBLIC anndatapp::anndatapp)
target_compile_options(ccell_preprocess PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS ccell_preprocess EXPORT ccellTargets
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_DATADIR}/ccell/preprocess/modules
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(CCELL_BUILD_TESTS)
  foreach(suite IN ITEMS numerics hvg dispatch transform)
    add_executable(ccell_preprocess_${suite} ${PROJECT_SOURCE_DIR}/tests/preprocess/test_${suite}.cpp)
    target_link_libraries(ccell_preprocess_${suite} PRIVATE ccell::preprocess)
    add_test(NAME ccell-preprocess-${suite} COMMAND ccell_preprocess_${suite})
  endforeach()
  add_executable(ccell_preprocess_zeroalloc ${PROJECT_SOURCE_DIR}/tests/preprocess/test_zeroalloc.cpp)
  target_link_libraries(ccell_preprocess_zeroalloc PRIVATE ccell::preprocess)
  if(CCELL_ENABLE_SANITIZERS)
    add_test(NAME ccell-preprocess-sanitized-behavior COMMAND ccell_preprocess_zeroalloc)
    set_tests_properties(ccell-preprocess-sanitized-behavior PROPERTIES SKIP_RETURN_CODE 77)
  else()
    add_test(NAME ccell-preprocess-zeroalloc COMMAND ccell_preprocess_zeroalloc)
    set_tests_properties(ccell-preprocess-zeroalloc PROPERTIES SKIP_RETURN_CODE 77)
  endif()

  add_executable(ccell_preprocess_example ${PROJECT_SOURCE_DIR}/examples/preprocess/basic.cpp)
  target_link_libraries(ccell_preprocess_example PRIVATE ccell::preprocess)
  add_test(NAME ccell-preprocess-example COMMAND ccell_preprocess_example)
endif()

if(CCELL_BUILD_QUALIFICATION)
  add_executable(ccell_preprocess_interop_driver ${PROJECT_SOURCE_DIR}/tests/preprocess/interop_driver.cpp)
  target_link_libraries(ccell_preprocess_interop_driver PRIVATE ccell::preprocess)
  add_executable(ccell_preprocess_benchmark ${PROJECT_SOURCE_DIR}/benchmarks/preprocess/qualification.cpp)
  target_link_libraries(ccell_preprocess_benchmark PRIVATE ccell::preprocess)
endif()
