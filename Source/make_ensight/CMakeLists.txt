
add_executable(make_ensight
  make_ensight.c
)
add_compile_options(-lm)
#target_link_libraries(<math.h> m)
target_link_libraries(make_ensight m)
