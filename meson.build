project('LSQuantumED', 'cpp', default_options : ['cpp_std=c++17', 'libdir=lib'], license : 'GPL-3.0', version : '0.0.1')

add_project_arguments(language: 'cpp')

conf_data = configuration_data()
conf_data.set_quoted('version', meson.project_version())
conf_data.set('DEVEL', get_option('DEVEL'))
conf_data.set('TIME_CODE', get_option('TIME_CODE'))
conf_data.set('BUILD_TEST', get_option('BUILD_TEST'))
conf_data.set('boostdir', get_option('boostdir'))
configure_file(input: 'include/meson_config.h.in',
	       output: 'meson_config.h',
	       install: true,
	       install_dir: 'include',
	       configuration: conf_data)

inc = include_directories(['.','include','include/spinOperators','include/tools'])

slepc_dep = dependency('SLEPc')

# Comment when running on marconi
mpi_dep = dependency('mpich')

subdir('include')
subdir('src')

subdir('examples')
subdir('postprocessing')

if get_option('BUILD_TEST') == true
   subdir('tests')
endif