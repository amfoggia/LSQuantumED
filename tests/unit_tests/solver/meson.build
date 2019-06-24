exec_solver_init = executable('solver_init.x',
                             ['test_solver_init.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                             include_directories: inc,
                             dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_solver1 = executable('solver1.x',
                         ['test_solver1.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_solver2 = executable('solver2.x',
                         ['test_solver2.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_solver3 = executable('solver3.x',
                         ['test_solver3.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_solver4 = executable('solver4.x',
                         ['test_solver4.cpp', env, time, rand_dis, lattice, basis, hamilt, solver],
                         include_directories: inc,
                         dependencies: [gtest_dep, petsc_dep, slepc_dep])

_tests = ['solver_init', 'solver1', 'solver2'] 

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
    tmp = get_variable('exec_' + _t)
    test(_t,
      mpiexec,
      args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t])
endforeach

_separate_tests = ['solver3', 'solver4']
foreach _t: _separate_tests
        tmp = get_variable('exec_' + _t)
    	test(_t,
     	mpiexec,
    	args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      	workdir: meson.current_source_dir(),
      	suite: ['unit_tests', _t])
endforeach
