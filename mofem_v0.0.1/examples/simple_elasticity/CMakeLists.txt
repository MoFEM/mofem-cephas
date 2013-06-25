# Lukasz Kaczmarczyk
# likask@civil.gla.ac.uk

add_executable(elasticity elasticity.cpp)
target_link_libraries(elasticity mofem_Field ${PROJECT_LIBS})

#example how to add new executable
#substitute "Program_Name" by name of new program
#add_executable(Program_Name Program_Name.cpp)
#target_link_libraries(Program_Name
#  fem ${PROJECT_LIBS} ${PROJECT_LIBS})

#add_executable(test_all test_all.cpp)
#target_link_libraries(test_all cgal_interface fem  ${PROJECT_LIBS} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
#target_link_libraries(test_all fem  ${PROJECT_LIBS})




