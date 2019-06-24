exec_env1 = executable('env1.x',
                      ['test_env1.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env2 = executable('env2.x',
                      ['test_env2.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env3 = executable('env3.x',
                      ['test_env3.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env4 = executable('env4.x',
                      ['test_env4.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env5 = executable('env5.x',
                      ['test_env5.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env6 = executable('env6.x',
                      ['test_env6.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env7 = executable('env7.x',
                      ['test_env7.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env8 = executable('env8.x',
                      ['test_env8.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env9 = executable('env9.x',
                      ['test_env9.cpp', env, time],
                      include_directories: inc,
                      dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env10 = executable('env10.x',
                       ['test_env10.cpp', env, time],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env11 = executable('env11.x',
                       ['test_env11.cpp', env, time],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_env12 = executable('env12.x',
                       ['test_env12.cpp', env, time],
                       include_directories: inc,
                       dependencies: [gtest_dep, petsc_dep, slepc_dep])

_tests = ['env'] 
_cases = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']

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
