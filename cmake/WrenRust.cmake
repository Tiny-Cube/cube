find_program(CARGO_EXECUTABLE cargo REQUIRED)
set(cargo_manifest "${CMAKE_CURRENT_SOURCE_DIR}/rust/Cargo.toml")
set(cargo_dir "${CMAKE_CURRENT_BINARY_DIR}/cargo")
set(cargo_profile "$<IF:$<CONFIG:Debug>,dev,release>")
set(cargo_output "${cargo_dir}/$<IF:$<CONFIG:Debug>,debug,release>")
set(cargo_flags --manifest-path "${cargo_manifest}" --target-dir "${cargo_dir}" --locked --offline)
file(GLOB_RECURSE rust_sources CONFIGURE_DEPENDS "rust/src/*.rs" "rust/examples/*.rs")
set(rust_dependencies "${cargo_manifest}" rust/Cargo.lock ${rust_sources})
set(rust_cli "${wren_bin}/ukmc-rust${CMAKE_EXECUTABLE_SUFFIX}")
set(rust_example "${wren_bin}/generic-rust${CMAKE_EXECUTABLE_SUFFIX}")
add_custom_command(OUTPUT "${rust_cli}" "${rust_example}"
  COMMAND "${CARGO_EXECUTABLE}" build ${cargo_flags} --profile "${cargo_profile}" --bin ukmc-rust --example generic
  COMMAND ${CMAKE_COMMAND} -E make_directory "${wren_bin}"
  COMMAND ${CMAKE_COMMAND} -E copy "${cargo_output}/ukmc-rust${CMAKE_EXECUTABLE_SUFFIX}" "${rust_cli}"
  COMMAND ${CMAKE_COMMAND} -E copy "${cargo_output}/examples/generic${CMAKE_EXECUTABLE_SUFFIX}" "${rust_example}"
  DEPENDS ${rust_dependencies} VERBATIM)
add_custom_target(wren_rust ALL DEPENDS "${rust_cli}" "${rust_example}")
install(PROGRAMS "${rust_cli}" DESTINATION "${CMAKE_INSTALL_BINDIR}")
set(rust_kernel "${wren_bin}/kernel-rust${CMAKE_EXECUTABLE_SUFFIX}")
add_custom_command(OUTPUT "${rust_kernel}"
  COMMAND "${CARGO_EXECUTABLE}" build ${cargo_flags} --profile "${cargo_profile}" --example kernel-rust
  COMMAND ${CMAKE_COMMAND} -E make_directory "${wren_bin}"
  COMMAND ${CMAKE_COMMAND} -E copy "${cargo_output}/examples/kernel-rust${CMAKE_EXECUTABLE_SUFFIX}" "${rust_kernel}"
  DEPENDS ${rust_dependencies} benchmarks/kernel.rs VERBATIM)
add_custom_target(wren_rust_benchmark DEPENDS "${rust_kernel}")
