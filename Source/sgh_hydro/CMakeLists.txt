
add_executable( sgh_hydro_3d
  driver.cpp
  )
add_compile_options(-O3 -lm -Wno-unused-result)
#target_link_libraries(<math.h> m)
target_link_libraries(sgh_hydro_3d m ${MPI_CXX_LIBRARIES} )
