exec_tools_random1 = executable('tools_random1.x',
                                ['test_tools_random1.cpp', rand_dis],
                                include_directories: inc,
                                dependencies: [gtest_dep, petsc_dep, slepc_dep])


_tests = ['tools_random'] 
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
