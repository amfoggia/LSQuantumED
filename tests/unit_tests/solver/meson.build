exec_solver_init = executable('solver_init.x',
                             ['test_solver_init.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                             include_directories: inc,
                             cpp_args:['-I'+boost],
			     dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_solver1 = executable('solver1.x',
                         ['test_solver1.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         cpp_args:['-I'+boost],
			 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_solver2 = executable('solver2.x',
                         ['test_solver2.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         cpp_args:['-I'+boost],
			 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_solver3 = executable('solver3.x',
                         ['test_solver3.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         cpp_args:['-I'+boost],
			 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_solver4 = executable('solver4.x',
                         ['test_solver4.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         cpp_args:['-I'+boost],
			 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_solver5 = executable('solver5.x',
                         ['test_solver5.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         cpp_args:['-I'+boost],
			 dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['_init', '1', '2', '5'] 

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
    tmp = get_variable('exec_solver' + _t)
    test('solver'+_t,
      mpiexec,
      args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', 'solver'])
endforeach

_separate_tests = ['3', '4']
foreach _t: _separate_tests
        tmp = get_variable('exec_solver' + _t)
    	test('solver'+_t,
     	mpiexec,
    	args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      	workdir: meson.current_source_dir(),
      	suite: ['unit_tests', 'solver'])
endforeach