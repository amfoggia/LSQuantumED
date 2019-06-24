exec_prealloc_info1 = executable('prealloc_info1.x',
                                 ['test_prealloc_info1.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_prealloc_info2 = executable('prealloc_info2.x',
                                 ['test_prealloc_info2.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_prealloc_info3 = executable('prealloc_info3.x',
                                 ['test_prealloc_info3.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

_tests = ['prealloc_info'] 
_cases = ['1', '3']

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

tmp = get_variable('exec_prealloc_info2')
test('prealloc_info',
     mpiexec,
     args: ['-np', '5', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'prealloc_info2'])