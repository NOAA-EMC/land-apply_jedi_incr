#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

#include "oops/base/Variables.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"

#include "noahmp_snow_increments.hpp"


namespace landincr {
  /**
   * AddSnowIncrement Class Implementation
   *
   */

class AddSnowIncrement : public oops::Application {
 public:
  using oops::Application::Application;
  explicit AddSnowIncrement(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}

  virtual ~AddLandIncrement() {}

  static const std::string classname() {
    return "landincr::AddSnowIncrement";
  }

  int execute(const eckit::Configuration& fullConfig) const override {
    // The actual State/Increment construction follows the existing
    // fv3jedi_add_soil_increment.h implementation.

    int myrank = this->getComm().rank();

    // We assume both state and increment are in the same geometry
    // const ijedi::Geometry<fv3jedi::Traits> geom_(eckit::LocalConfiguration(fullConfig, "geometry"),this->getComm());
    fv3jedi::Geometry geom_(eckit::LocalConfiguration(fullConfig, "geometry"),this->getComm());
    // oops::Log::info() << "geom ny "<<geom_.npy() << " nx " << geom_.npx() << std::endl;

    // Read state
    // ijedi::State<ijedi::Traits> xx(geom_, eckit::LocalConfiguration(fullConfig, "background state"));
    fv3jedi::State xx(geom_, eckit::LocalConfiguration(fullConfig, "background state"));
    oops::Log::test() << "Background state: " << xx << std::endl;

    // Read increment
    // oops::Log::info() << "Reading increment" << std::endl;
    const eckit::LocalConfiguration incParams(fullConfig, "increment");
    oops::Variables incVars(incParams, "variables");
    // ijedi::Increment<ijedi::Traits> dx(geom_, incVars, xx.validTime());
    fv3jedi::Increment dx(geom_, incVars, xx.validTime());
    dx.read(incParams);
    oops::Log::test() << "Increment: " << dx << std::endl;

    atlas::FieldSet bkg_fs;
    atlas::FieldSet inc_fs;

    xx.toFieldSet(bkg_fs);
    dx.toFieldSet(inc_fs);

    if (!bkg_fs.has("sheleg") ||
        !bkg_fs.has("snwdph") ||
        !bkg_fs.has("snowxy") ||
        !bkg_fs.has("sneqvoxy") ||
        !bkg_fs.has("zsnsoxy") ||
        !bkg_fs.has("tsnoxy") ||
        !bkg_fs.has("snicexy") ||
        !bkg_fs.has("snliqxy") ||
        !bkg_fs.has("stc")) {
      throw eckit::BadValue(
          "Required Noah-MP snow fields are missing", Here());
    }

    auto swe_view =
        atlas::array::make_view<double, 2>(bkg_fs["sheleg"]);

    auto snow_depth_view =
        atlas::array::make_view<double, 2>(bkg_fs["snwdph"]);

    auto active_layers_view =
        atlas::array::make_view<double, 2>(bkg_fs["snowxy"]);

    auto swe_previous_view =
        atlas::array::make_view<double, 2>(bkg_fs["sneqvoxy"]);

    auto interface_view =
        atlas::array::make_view<double, 3>(bkg_fs["zsnsoxy"]);

    auto snow_temperature_view =
        atlas::array::make_view<double, 3>(bkg_fs["tsnoxy"]);

    auto snow_ice_view =
        atlas::array::make_view<double, 3>(bkg_fs["snicexy"]);

    auto snow_liq_view =
        atlas::array::make_view<double, 3>(bkg_fs["snliqxy"]);

    auto soil_temperature_view =
        atlas::array::make_view<double, 2>(bkg_fs["stc"]);

    const std::size_t nland = swe_view.shape(0);

    snow_increments::NoahMPSnowState snow_state;
    snow_state.resize(nland);

    for (std::size_t i = 0; i < nland; ++i) {
      snow_state.swe[i] = swe_view(i, 0);
      snow_state.snow_depth[i] = snow_depth_view(i, 0);
      snow_state.active_snow_layers[i] = active_layers_view(i, 0);
      snow_state.swe_previous[i] = swe_previous_view(i, 0);

      // Noah-MP snow arrays have three active snow layers.
      for (std::size_t k = 0; k < snow_increments::kSnowLayers; ++k) {
        snow_state.temperature_snow[i][k] = snow_temperature_view(i, k);

        snow_state.snow_ice_layer[i][k] = snow_ice_view(i, k);

        snow_state.snow_liq_layer[i][k] = snow_liq_view(i, k);
      }

      // zsnsoxy has seven interfaces.
      for (std::size_t k = 0; k < snow_increments::kSnowInterfaces; ++k) {
        snow_state.snow_soil_interface[i][k] = interface_view(i, k);
      }

      // The Fortran code uses the first soil-temperature layer.
      snow_state.temperature_soil[i] = soil_temperature_view(i, 0);
    }

    std::vector<double> snow_depth_increment(nland, 0.0);

    if (!inc_fs.has("snwdph")) {
      throw eckit::BadValue("Snow-depth increment field snwdph is missing", Here());
    }

    auto increment_view = atlas::array::make_view<double, 2>(inc_fs["snwdph"]);

    for (std::size_t i = 0; i < nland; ++i) {
      snow_depth_increment[i] = increment_view(i, 0);
    }

    double noincr_threshold = 999999999.9;
    bool print_summary = false;
    bool print_debug = false;

    incParams.get("noincr_threshold", noincr_threshold);
    incParams.get("print_summary", print_summary);
    incParams.get("print_debug", print_debug);

    snow_increments::NoahMPSnowIncrementer::updateAllLayers(
        nland,
        snow_depth_increment,
        snow_state,
        noincr_threshold,
        print_summary,
        print_debug);

    // Copy the updated land-vector state back to Atlas.
    for (std::size_t i = 0; i < nland; ++i) {
      swe_view(i, 0) = snow_state.swe[i];
      snow_depth_view(i, 0) = snow_state.snow_depth[i];
      active_layers_view(i, 0) = snow_state.active_snow_layers[i];
      swe_previous_view(i, 0) = snow_state.swe_previous[i];

      for (std::size_t k = 0; k < snow_increments::kSnowLayers; ++k) {
        snow_temperature_view(i, k) = snow_state.temperature_snow[i][k];

        snow_ice_view(i, k) = snow_state.snow_ice_layer[i][k];

        snow_liq_view(i, k) = snow_state.snow_liq_layer[i][k];
      }

      for (std::size_t k = 0; k < snow_increments::kSnowInterfaces; ++k) {
        interface_view(i, k) = snow_state.snow_soil_interface[i][k];
      }
    }


    xx.fromFieldSet(bkg_fs);
    oops::Log::test() << "Updated State: " << xx << std::endl;

    // Write updated state to file
    xx.write(eckit::LocalConfiguration(fullConfig, "output state"));

    return 0;
  }

 private:
   // hard coded defaults--unlikely to change
   static constexpr int veg_type_landice = 15;
   bool frac_grid = true;
   float fice_threshold = 0.0;
   float lfrac_threshold = 0.0001;

   std::string appname() const override {
     return "landincr::AddSnowIncrement";
   }
};

}  // namespace landincr
