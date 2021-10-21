find_package(precice REQUIRED CONFIG)

set(programs
		dune-perpendicular-flap
    dune-perpendicular-flap-implicit)

#include precice
find_package(precice)

foreach(_program ${programs})
	add_executable(${_program} ${_program}.cc)
	target_link_libraries(${_program} PRIVATE precice::precice)
endforeach()
