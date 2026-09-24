################################################################################
# ijedi_soil_incr
################################################################################

list( APPEND sorc_files_ijedisoilinc
    ./ijedi_add_soil_increment.cc
    ./ijedi_add_soil_increment.h
    ./soil_increments.hpp
)

ecbuild_add_executable( TARGET          ijedi_soil_incr.x
                        SOURCES         ${sorc_files_ijedisoilinc} )
target_compile_features( ijedi_soil_incr.x PUBLIC cxx_std_17)
target_link_libraries( ijedi_soil_incr.x PUBLIC NetCDF::NetCDF_CXX NetCDF::NetCDF_Fortran oops atlas ijedi)

################################################################################
