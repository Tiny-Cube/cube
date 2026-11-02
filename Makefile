.PHONY: all native test package
all: native
native:
	python3 scripts/build.py
test:
	python3 scripts/check.py
package:
	uv build
