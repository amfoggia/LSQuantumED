exec_lattice1 = executable('lattice1.x',
                          ['test_lattice1.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice2 = executable('lattice2.x',
                          ['test_lattice2.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice3 = executable('lattice3.x',
                          ['test_lattice3.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice4 = executable('lattice4.x',
                          ['test_lattice4.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice5 = executable('lattice5.x',
                          ['test_lattice5.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice6 = executable('lattice6.x',
                          ['test_lattice6.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice7 = executable('lattice7.x',
                          ['test_lattice7.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])			  

exec_lattice8 = executable('lattice8.x',
                          ['test_lattice8.cpp', env, time, lattice],
                          include_directories: inc,
			  cpp_args:['-I'+boost],
                          dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_lattice10 = executable('lattice10.x',
                           ['test_lattice10.cpp', env, time, lattice],
                           include_directories: inc,
		           cpp_args:['-I'+boost],
                           dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_lattice11 = executable('lattice11.x',
                           ['test_lattice11.cpp', env, time, lattice],
                           include_directories: inc,
		           cpp_args:['-I'+boost],
                           dependencies: [gtest_dep, slepc_dep, mpi_dep])

# exec_lattice9 = executable('lattice9.x',
#                           ['test_lattice9.cpp', env, time, lattice],
#                           include_directories: inc,
#                           cpp_args:['-I'+boost],
#                           dependencies: [gtest_dep, slepc_dep, mpi_dep])

_tests = ['lattice'] 
_cases = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11']

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
