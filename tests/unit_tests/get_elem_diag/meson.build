exec_get_elem_diag1 = executable('get_elem_diag1.x',
                                ['test_get_elem_diag1.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
				cpp_args:['-I'+boost],
                                dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_elem_diag2 = executable('get_elem_diag2.x',
                                ['test_get_elem_diag2.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
				cpp_args:['-I'+boost],
                                dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_elem_diag3 = executable('get_elem_diag3.x',
                                ['test_get_elem_diag3.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
				cpp_args:['-I'+boost],
                                dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_elem_diag4 = executable('get_elem_diag4.x',
                                ['test_get_elem_diag4.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
				cpp_args:['-I'+boost],
                                dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_elem_diag5 = executable('get_elem_diag5.x',
                                ['test_get_elem_diag5.cpp', env, time, lattice, basis, hamilt],
                                include_directories: inc,
				cpp_args:['-I'+boost],
                                dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['get_elem_diag'] 
_cases = ['1', '2', '3', '4', '5']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  foreach _num: _cases
    tmp = get_variable('exec_' + _t + _num)
    test(_t+_num,
      mpiexec,
      args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t])
    endforeach
endforeach
