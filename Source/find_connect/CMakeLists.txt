
add_executable( find_connect
  find_connect.cpp
)
add_compile_options(-lm)
#target_link_libraries(<math.h> m)
target_link_libraries(find_connect m)
