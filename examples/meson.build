install_subdir('.',
	       install_dir: 'examples',
	       strip_directory: true,
	       exclude_files: ['meson.build'])

examples_src = [env,
               basis,
	       lattice,
	       sublat,
	       hamilt,
	       solver,
	       rand_dis,
	       time,
	       correl]

exec_main_NO_dis = executable('main_NO_disorder.x',
                             ['main_NO_disorder.cpp', examples_src],
                             include_directories: inc,
			     cpp_args: ['-I'+boost],
                             dependencies: [slepc_dep, mpi_dep],
			     install: true,
                             install_dir: 'examples')

exec_main_dis = executable('main_disorder.x',
                          ['main_disorder.cpp', examples_src],
                          include_directories: inc,
			  cpp_args: ['-I'+boost],
                          dependencies: [slepc_dep, mpi_dep],
                          install: true,
		          install_dir: 'examples')
			  

exec_main = executable('main.x',
                      ['main.cpp', examples_src],
                      include_directories: inc,
		      cpp_args: ['-I'+boost],
                      dependencies: [slepc_dep, mpi_dep],
		      install: true,
		      install_dir: 'examples')

exec_ex1 = executable('ex1.x',
                     ['ex1.cpp', examples_src],
                     include_directories: inc,
	             cpp_args: ['-I'+boost],
                     dependencies: [slepc_dep, mpi_dep],
	             install: true,
	             install_dir: 'examples')

exec_ex2 = executable('ex2.x',
                     ['ex2.cpp', examples_src],
                     include_directories: inc,
	             cpp_args: ['-I'+boost],
                     dependencies: [slepc_dep, mpi_dep],
	             install: true,
	             install_dir: 'examples')

exec_ex3 = executable('ex3.x',
                     ['ex3.cpp', examples_src],
                     include_directories: inc,
	             cpp_args: ['-I'+boost],
                     dependencies: [slepc_dep, mpi_dep],
	             install: true,
	             install_dir: 'examples')

exec_ex4 = executable('ex4.x',
                     ['ex4.cpp', examples_src],
                     include_directories: inc,
	             cpp_args: ['-I'+boost],
                     dependencies: [slepc_dep, mpi_dep],
	             install: true,
	             install_dir: 'examples')

mpiexec = find_program('mpiexec')

run_target('example1',
	  command: [mpiexec, '-np', '1', exec_main, '-n 6'])

run_target('example2',
	  command: [mpiexec, '-np', '1', exec_main_dis, '-n 6'])

run_target('example3',
	  command: [mpiexec, '-np', '1', exec_main_NO_dis, '-n 6'])