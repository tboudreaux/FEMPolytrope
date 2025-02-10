coeff_sources = files(
  'coeff.cpp'
)

coeff_headers = files(
  'coeff.h'
)

libcoeff = static_library('coeff',
  coeff_sources,
  include_directories : include_directories('.'),
  cpp_args: ['-fvisibility=default'],
  dependencies: [mfem_dep],
  install: true
)

install_headers(coeff_headers, subdir : 'coeff')
