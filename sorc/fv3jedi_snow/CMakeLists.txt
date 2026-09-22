################################################################################
# fv3jedi_snow_incr
################################################################################

list(APPEND sorc_files_snowinc
    ./fv3jedi_add_snow_increment.cc
    ./fv3jedi_add_snow_increment.h
    ./noahmp_snow_increments.hpp
)

ecbuild_add_executable(
    TARGET fv3jedi_snow_incr.x
    SOURCES ${sorc_files_snowinc}
)

target_compile_features(
    fv3jedi_snow_incr.x PUBLIC cxx_std_17
)

target_link_libraries(
    fv3jedi_snow_incr.x PUBLIC
    NetCDF::NetCDF_CXX
    NetCDF::NetCDF_Fortran
    oops
    atlas
    fv3jedi
)

################################################################################
