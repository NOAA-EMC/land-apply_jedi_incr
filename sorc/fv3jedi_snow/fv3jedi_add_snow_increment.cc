#include "fv3jedi_add_snow_increment.h"

#include "fv3jedi/Utilities/Traits.h"

#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  //AddLandIncrement<ijedi::Traits> addLandIncrement;
  landincr::AddSnowIncrement addSnowIncrement;
  return run.execute(addSnowIncrement);
}
