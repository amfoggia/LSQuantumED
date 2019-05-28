exec_szsz1 = executable('szsz1.x',
                       ['test_szsz1.cpp', env, time, lattice, basis, hamilt, solver, correl],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

exec_szsz2 = executable('szsz2.x',
                       ['test_szsz2.cpp', env, time, lattice, basis, hamilt, solver, correl],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

# exec_szsz3 = executable('szsz3.x',
#                                 ['test_szsz3.cpp', env, time, lattice, basis, hamilt],
#                                 include_directories: inc,
#                                 dependencies: [gtest_dep, petsc_dep, slepc_dep, mpi_dep])

_tests = ['szsz'] 
_cases = ['1']

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

_tests = ['szsz'] 
_cases = ['2']

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

foreach _t: _tests
  foreach _num: _cases
    tmp = get_variable('exec_' + _t + _num)
    test(_t,
      mpiexec,
      args: ['-np', '4', tmp.full_path(), '-malloc_dump'],
      workdir: meson.current_source_dir(),
      suite: ['unit_tests', _t+_num])
    endforeach
endforeach
