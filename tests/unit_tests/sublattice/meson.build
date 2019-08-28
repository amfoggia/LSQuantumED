exec_sublattice1 = executable('sublattice1.x',
                             ['test_sublattice1.cpp', env, time, lattice, sublat],
                             include_directories: inc,
			     cpp_args:['-I'+boost],
                             dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_sublattice2 = executable('sublattice2.x',
                             ['test_sublattice2.cpp', env, time, lattice, sublat],
                             include_directories: inc,
			     cpp_args:['-I'+boost],
                             dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_sublattice3 = executable('sublattice3.x',
                             ['test_sublattice3.cpp', env, time, lattice, sublat],
                             include_directories: inc,
			     cpp_args:['-I'+boost],
                             dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['sublattice'] 
_cases = ['1','2', '3']

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
