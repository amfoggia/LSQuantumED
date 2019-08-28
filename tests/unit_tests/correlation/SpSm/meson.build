exec_spsm1 = executable('spsm1.x',
                       ['test_spsm1.cpp', env, time, lattice, basis, hamilt, solver, correl],
                       include_directories: inc,
		       cpp_args:['-I'+boost],
                       dependencies: [gtest_dep, slepc_dep, mpi_dep])

# exec_spsm2 = executable('spsm2.x',
#                                 ['test_spsm2.cpp', env, time, lattice, basis, hamilt],
#                                 include_directories: inc,
#                                 cpp_args:['-I'+boost],
#                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

# exec_spsm3 = executable('spsm3.x',
#                                 ['test_spsm3.cpp', env, time, lattice, basis, hamilt],
#                                 include_directories: inc,
#                                 cpp_args:['-I'+boost],
#                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['spsm'] 
_cases = ['1']

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