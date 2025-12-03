project(
    'highfive', 'cpp',
    version : 'v2.10.1',
)

cxx = meson.get_compiler('cpp')

includes = include_directories('include', is_system: not meson.is_subproject())

# find HDF5 (or serial HDF5 if any) & add it to dependency
project_deps = []
h5dep = dependency('hdf5', language: 'cpp', required: false)
if h5dep.found()
    project_deps += h5dep
else
    h5lib = cxx.find_library('hdf5', required: false)
    if not h5lib.found()
        project_deps += cxx.find_library('hdf5_serial', required: true)
    else
        project_deps += h5lib
    endif
endif

highfive_dep = declare_dependency(
    include_directories: includes,
    dependencies: project_deps
)