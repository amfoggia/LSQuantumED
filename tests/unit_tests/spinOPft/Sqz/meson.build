exec_sqz1 = executable('sqz1.x',
                       ['test_sqz1.cpp', env, time, lattice, basis, hamilt, solver],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

# exec_sqz2 = executable('sqz2.x',
#                        ['test_sqz2.cpp', env, time, lattice, basis, hamilt],
#                        include_directories: inc,
#                        dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

_tests = ['sqz'] 
_cases = ['1']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  foreach _num: _cases
    tmp = get_variable('exec_' + _t + _num)
    test(_t,
      mpiexec,
      args: ['-np', '3', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t+_num])
    endforeach
endforeach
