find_program(UV_EXECUTABLE uv REQUIRED)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/pyproject.toml")
execute_process(COMMAND "${Python3_EXECUTABLE}" -c
  "import sys,tomllib; p=tomllib.load(open(sys.argv[1],'rb'))['project']; print(';'.join([p['version'],*p['scripts']]))"
  "${PROJECT_SOURCE_DIR}/pyproject.toml" OUTPUT_VARIABLE wren_entry_points
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
list(POP_FRONT wren_entry_points wren_python_version)
if(NOT wren_python_version STREQUAL PROJECT_VERSION)
  message(FATAL_ERROR "Wren CMake version ${PROJECT_VERSION} differs from pyproject.toml ${wren_python_version}")
endif()
set(wren_dist "${CMAKE_CURRENT_BINARY_DIR}/dist")
set(wren_wheel "${wren_dist}/wren-${PROJECT_VERSION}-py3-none-any.whl")
set(wren_stage "${CMAKE_CURRENT_BINARY_DIR}/python-package")
file(GLOB_RECURSE wren_python_sources CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/wren/*.py")
add_custom_command(OUTPUT "${wren_wheel}" "${wren_stage}/wren/__init__.py"
  COMMAND "${UV_EXECUTABLE}" build --wheel --out-dir "${wren_dist}" "${PROJECT_SOURCE_DIR}"
  COMMAND ${CMAKE_COMMAND} -E rm -rf "${wren_stage}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${wren_stage}"
  COMMAND "${Python3_EXECUTABLE}" -m zipfile -e "${wren_wheel}" "${wren_stage}"
  DEPENDS pyproject.toml README.md LICENSE uv.lock ${wren_python_sources}
  COMMENT "Building the Wren Python wheel with uv" VERBATIM)
add_custom_target(wren_python ALL DEPENDS "${wren_wheel}" "${wren_stage}/wren/__init__.py")
add_custom_target(wren_sdist
  COMMAND "${UV_EXECUTABLE}" build --sdist --out-dir "${wren_dist}" "${PROJECT_SOURCE_DIR}"
  COMMENT "Building the Wren source distribution with uv" VERBATIM)
add_custom_target(wren_package DEPENDS wren_python wren_sdist)

# Install the wheel for uv users and a relocatable unpacked copy for the CLI.
# Unpacked files are needed for the runtime's source/provenance hashes.
install(FILES "${wren_wheel}" DESTINATION "${CMAKE_INSTALL_DATADIR}/wren/wheels")
install(DIRECTORY "${wren_stage}/" DESTINATION "${CMAKE_INSTALL_DATADIR}/wren/python")
file(RELATIVE_PATH WREN_PYTHON_FROM_BIN "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/wren/python")
configure_file("${PROJECT_SOURCE_DIR}/cmake/wren-launcher.py.in" "${CMAKE_CURRENT_BINARY_DIR}/wren-launcher" @ONLY)
foreach(entry IN LISTS wren_entry_points)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/wren-launcher" DESTINATION "${CMAKE_INSTALL_BINDIR}" RENAME "${entry}")
endforeach()
