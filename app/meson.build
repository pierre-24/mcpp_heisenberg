app_files = files(
    'app.cpp',
    'main.cpp',
)

# add dependencies
app_deps = [lmcpp_heisenberg_dep]

includes = include_directories('.')

cli11_proj = subproject('cli11')
app_deps += cli11_proj.get_variable('CLI11_dep')

tomlplusplus_proj = subproject('tomlplusplus')
app_deps += tomlplusplus_proj.get_variable('tomlplusplus_dep')

# make executable
mcpp_heisenberg_app = executable(
    meson.project_name() + '_run',
    app_files,
    include_directories : includes,
    dependencies: app_deps,
    cpp_args : build_args,
    install : true
)