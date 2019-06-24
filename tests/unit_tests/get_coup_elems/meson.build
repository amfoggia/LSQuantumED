exec_get_coup_elems1 = executable('get_coup_elems1.x',
                                 ['test_get_coup_elems1.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_get_coup_elems2 = executable('get_coup_elems2.x',
                                 ['test_get_coup_elems2.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

exec_get_coup_elems3 = executable('get_coup_elems3.x',
                                 ['test_get_coup_elems3.cpp', env, time, lattice, basis, hamilt],
                                 include_directories: inc,
                                 dependencies: [gtest_dep, petsc_dep, slepc_dep])

_prefix = meson.current_source_dir()

mpiexec = find_program('mpiexec')

tmp = get_variable('exec_' + 'get_coup_elems' + '1')
test('get_coup_elems',
     mpiexec,
     args: ['-np', '2', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems1'])

tmp = get_variable('exec_' + 'get_coup_elems' + '2')
test('get_coup_elems',
     mpiexec,
     args: ['-np', '3', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
     suite: ['unit_tests', 'get_coup_elems2'])

tmp = get_variable('exec_' + 'get_coup_elems' + '3')
test('get_coup_elems',
     mpiexec,
     args: ['-np', '5', tmp.full_path(), '-malloc_dump'],
     workdir: meson.current_source_dir(),
suite: ['unit_tests', 'get_coup_elems3'])