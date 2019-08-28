exec_sqz1 = executable('sqz1.x',
                       ['test_sqz1.cpp', env, time, lattice, basis, hamilt, solver],
                       include_directories: inc,
                       cpp_args: ['-I'+boost],
                       dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_sqp1 = executable('sqp1.x',
                       ['test_sqp1.cpp', env, time, lattice, basis, hamilt],
                       include_directories: inc,
		       cpp_args: ['-I'+boost],
                       dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['sq'] 
_cases = ['z1','p1']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  foreach _num: _cases
    tmp = get_variable('exec_' + _t + _num)
    test(_t+_num,
      mpiexec,
      args: ['-np', '3', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t])
    endforeach
endforeach
