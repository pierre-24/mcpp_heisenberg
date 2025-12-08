app_lib_files = files(
    'app.cpp',
)

# add dependencies
app_deps = [lmcpp_heisenberg_dep]

includes = include_directories('.')

cli11_proj = subproject('cli11')
app_deps += cli11_proj.get_variable('CLI11_dep')

tomlplusplus_proj = subproject('tomlplusplus')
app_deps += tomlplusplus_proj.get_variable('tomlplusplus_dep')

# make library for runner
mcpp_heisenberg_run_lib = library(
    'lib' + meson.project_name() + '_run',
    app_lib_files,
    install : true,
    cpp_args : build_args,
    include_directories : includes,
    dependencies: app_deps,
)

lmcpp_heisenberg_run_dep = declare_dependency(
    link_with: mcpp_heisenberg_run_lib,
    include_directories: includes,
    dependencies: app_deps
)

# make executable
mcpp_heisenberg_app = executable(
    meson.project_name() + '_run',
    'main.cpp',
    include_directories : includes,
    dependencies: lmcpp_heisenberg_run_dep,
    cpp_args : build_args,
    install : true
)