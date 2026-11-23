BUILD_DIR ?= build/cmake
.PHONY: all native test package
all: native
native:
	cmake -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE=Release
	cmake --build "$(BUILD_DIR)" --config Release --parallel
test: native
	ctest --test-dir "$(BUILD_DIR)" -C Release --output-on-failure
package:
	cmake -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE=Release
	cmake --build "$(BUILD_DIR)" --config Release --target wren_package
