exec_basis1 = executable('basis1.x',
                        ['test_basis1.cpp', env, time, basis],
                        include_directories: inc,
			cpp_args:['-I'+boost],
                        dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_basis2 = executable('basis2.x',
                        ['test_basis2.cpp', env, time, basis],
                        include_directories: inc,
			cpp_args:['-I'+boost],
                        dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_basis3 = executable('basis3.x',
                        ['test_basis3.cpp', env, time, basis],
                        include_directories: inc,
			cpp_args:['-I'+boost],
                        dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_basis4 = executable('basis4.x',
                        ['test_basis4.cpp', env, time, basis],
                        include_directories: inc,
			cpp_args:['-I'+boost],
                        dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['basis'] 
_cases = ['1', '2']

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

tmp = get_variable('exec_' + 'basis' + '3')
test('basis3',
     mpiexec,
     args: ['-np', '1', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'basis'])
     
tmp = get_variable('exec_' + 'basis' + '4')
test('basis4',
     mpiexec,
     args: ['-np', '5', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'basis'])
