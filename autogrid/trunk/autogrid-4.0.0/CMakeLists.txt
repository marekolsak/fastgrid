PROJECT(autogrid)

SET(AUTOGRID_SOURCES
main.cpp
inputdataloader.cpp
logfile.cpp
pairwiseinteractionenergies.cpp
times.cpp
programparameters.cpp
utils.cpp
bondvectors.cpp
parameterlibrary.cpp
desolvexpfunc.cpp
gridmap.cpp
)

ADD_EXECUTABLE(autogrid4 ${AUTOGRID_SOURCES})
