exec_get_coup_elems1 = executable('get_coup_elems1.x',
                                 ['test_get_coup_elems1.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
				 cpp_args:['-I'+boost],
                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_coup_elems2 = executable('get_coup_elems2.x',
                                 ['test_get_coup_elems2.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
				 cpp_args:['-I'+boost],
                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_coup_elems3 = executable('get_coup_elems3.x',
                                 ['test_get_coup_elems3.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
				 cpp_args:['-I'+boost],
                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_coup_elems4 = executable('get_coup_elems4.x',
                                 ['test_get_coup_elems4.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
				 cpp_args:['-I'+boost],
                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

exec_get_coup_elems5 = executable('get_coup_elems5.x',
                                 ['test_get_coup_elems5.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
				 cpp_args:['-I'+boost],
                                 dependencies: [gtest_dep, slepc_dep, mpi_dep])

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

tmp = get_variable('exec_' + 'get_coup_elems' + '1')
test('get_coup_elems1',
     mpiexec,
     args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems'])

tmp = get_variable('exec_' + 'get_coup_elems' + '2')
test('get_coup_elems2',
     mpiexec,
     args: ['-np', '3', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems'])

tmp = get_variable('exec_' + 'get_coup_elems' + '3')
test('get_coup_elems3',
     mpiexec,
     args: ['-np', '5', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
suite: ['unit_tests', 'get_coup_elems'])

tmp = get_variable('exec_' + 'get_coup_elems' + '4')
test('get_coup_elems4',
     mpiexec,
     args: ['-np', '1', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems'])

tmp = get_variable('exec_' + 'get_coup_elems' + '5')
test('get_coup_elems5',
     mpiexec,
     args: ['-np', '1', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems'])
