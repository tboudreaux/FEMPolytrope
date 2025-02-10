utils_sources = files(
  'io.cpp',
  'mfemUtils.cpp'
)

utils_headers = files(
  'io.h',
  'mfemUtils.h'
)

libutils = static_library('utils',
  utils_sources,
  include_directories : include_directories('.'),
  cpp_args: ['-fvisibility=default'],
  dependencies: [mfem_dep],
  install: true
)

install_headers(utils_headers, subdir : 'utils')
