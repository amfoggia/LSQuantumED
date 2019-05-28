exec_sublattice1 = executable('sublattice1.x',
                             ['test_sublattice1.cpp', env, time, lattice, sublat],
                             include_directories: inc,
                             dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

exec_sublattice2 = executable('sublattice2.x',
                             ['test_sublattice2.cpp', env, time, lattice, sublat],
                             include_directories: inc,
                             dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])
			  
_tests = ['sublattice'] 
_cases = ['1','2']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  foreach _num: _cases
    tmp = get_variable('exec_' + _t + _num)
    test(_t,
      mpiexec,
      args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t+_num])
    endforeach
endforeach
