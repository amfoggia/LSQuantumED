exec_get_disorder_elem1 = executable('get_disorder_elem1.x',
                                ['test_get_disorder_elem1.cpp', env, time, rand_dis, lattice, basis, hamilt],
                                include_directories: inc,
                                dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

exec_get_disorder_elem3 = executable('get_disorder_elem3.x',
                                ['test_get_disorder_elem3.cpp', env, time, rand_dis, lattice, basis, hamilt],
                                include_directories: inc,
                                dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

_tests = ['get_disorder_elem'] 
_cases = ['1','3']

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
