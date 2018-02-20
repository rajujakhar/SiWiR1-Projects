CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( heat )

SET( SRCS

   heattransfer.cpp
   ParabolicSolver.cpp
   Implicit.cpp
   Explicit.cpp
   global.cpp
   #mpi_cart_Create.cpp
   # more *.cpp files follow here
)


#SET(CMAKE_CXX_COMPILER mpicxx)

#add_compile_options(-std=c++0x) 

#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3" )

IF( CMAKE_COMPILER_IS_GNUCC )
   #IF( ENABLE_ANSI )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} —std=c++11 ansi -pedantic" )
   #ENDIF()
ENDIF()

ADD_EXECUTABLE( heat ${SRCS} )
