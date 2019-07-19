basis = files('basis.cpp')

lattice = files('lattice.cpp')

sublat = files('sublattice.cpp')

hamilt = files('hamiltonian.cpp')

solver = files('solver.cpp')

env = files('environment.cpp')

phys = files('phys.cpp')

subdir('spinOperators')

subdir('tools')

src_files = [basis,
	     lattice,
	     sublat,
	     hamilt,
	     solver,
	     env,
	     phys,
	     src_spinOperators,
	     src_tools]

boost = get_option('boostdir')

libLSQuantumED = build_target('lsquantumed',
			      src_files,
			      name_suffix: [],
			      cpp_args: ['-O2', '-Wall', '-I'+boost],
			      include_directories: inc,
			      dependencies: [slepc_dep, mpi_dep],
			      target_type: 'shared_library',
			      install: true)

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries: libLSQuantumED,
                 version: meson.project_version(),
		 name: 'LSQuantumED',
		 filebase: 'LSQuantumED',
		 extra_cflags: ['-O3', '-Wall', '-I'+boost],
		 requires: [slepc_dep, mpi_dep],
 		 description: 'Library to simulate spin systems on lattices with nearest and next-nearest neighbours interactions.')

thislib_dep = declare_dependency(dependencies: [slepc_dep, mpi_dep],
	      		         include_directories: inc,
				 link_args: ['-I'+boost],
				 link_with: [libLSQuantumED],
				 version: meson.project_version())