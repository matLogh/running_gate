# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})


add_library(rgfitters STATIC Fitters.h Fitters.cpp)
target_link_libraries(rgfitters PUBLIC theuerkauf_fitter ${ROOT_LIBRARIES})

# Add the RunningGate library
add_library(RunningGate STATIC
    RunningGate.cpp
    RunningGate.h
)
target_include_directories(RunningGate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(RunningGate PUBLIC rgfitters ${ROOT_LIBRARIES} )
