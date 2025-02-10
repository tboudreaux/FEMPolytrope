project('FEM1D', 'cpp')

cpp_compiler = meson.get_compiler('cpp')

mfem_inc = include_directories('/home/tboudreaux/.local/mfem/include')
mfem_lib = cpp_compiler.find_library('mfem', dirs: '/home/tboudreaux/.local/mfem/lib', required: true)
mfem_dep = declare_dependency(dependencies: mfem_lib, include_directories: mfem_inc)

subdir('src')
