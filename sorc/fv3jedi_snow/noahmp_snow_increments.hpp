#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace snow_increments {

// Noah-MP fixed dimensions.
constexpr std::size_t kSnowInterfaces = 7;
constexpr std::size_t kSnowLayers = 3;

// Interface values used by Noah-MP when no active snow layer exists.
constexpr std::array<double, kSnowInterfaces> kZeroSnowInterfaces = {
    0.0, 0.0, 0.0, -0.1, -0.4, -1.0, -2.0
};

// Soil-interface depths used when reconstructing the lower interfaces.
constexpr std::array<double, kSnowInterfaces> kSoilInterfaces = {
    0.0, 0.0, 0.0, 0.1, 0.4, 1.0, 2.0
};

/**
 * @brief Noah-MP snow state stored over land points only.
 *
 * The first dimension is the land-vector index, not the full FV3 tile index.
 *
 * snow_soil_interface:
 *   [land point][0..6]
 *
 * temperature_snow, snow_ice_layer, snow_liq_layer:
 *   [land point][0..2]
 */
struct NoahMPSnowState {
  std::vector<double> swe;
  std::vector<double> snow_depth;
  std::vector<double> active_snow_layers;
  std::vector<double> swe_previous;

  std::vector<std::array<double, kSnowInterfaces>> snow_soil_interface;
  std::vector<std::array<double, kSnowLayers>> temperature_snow;
  std::vector<std::array<double, kSnowLayers>> snow_ice_layer;
  std::vector<std::array<double, kSnowLayers>> snow_liq_layer;

  std::vector<double> temperature_soil;

  std::string name_snow_depth;
  std::string name_swe;

  /**
   * @brief Resize and initialize a Noah-MP snow state.
   */
  void resize(std::size_t npoints) {
    swe.assign(npoints, 0.0);
    snow_depth.assign(npoints, 0.0);
    active_snow_layers.assign(npoints, 0.0);
    swe_previous.assign(npoints, 0.0);

    snow_soil_interface.assign(
        npoints, kZeroSnowInterfaces);

    temperature_snow.assign(
        npoints, std::array<double, kSnowLayers>{});

    snow_ice_layer.assign(
        npoints, std::array<double, kSnowLayers>{});

    snow_liq_layer.assign(
        npoints, std::array<double, kSnowLayers>{});

    temperature_soil.assign(npoints, 0.0);
  }

  [[nodiscard]] std::size_t size() const {
    return swe.size();
  }

  [[nodiscard]] bool consistent() const {
    const std::size_t n = swe.size();

    return snow_depth.size() == n &&
           active_snow_layers.size() == n &&
           swe_previous.size() == n &&
           snow_soil_interface.size() == n &&
           temperature_snow.size() == n &&
           snow_ice_layer.size() == n &&
           snow_liq_layer.size() == n &&
           temperature_soil.size() == n;
  }
};

/**
 * @brief Summary statistics from one snow increment update.
 */
struct SnowIncrementSummary {
  std::size_t no_increment = 0;
  std::size_t removed_all_snow = 0;
  std::size_t added_multilayer = 0;
  std::size_t removed_multilayer = 0;
  std::size_t added_zero_layer = 0;
  std::size_t skipped_too_warm = 0;
  std::size_t created_layer = 0;
  std::size_t removed_zero_layer = 0;
  std::size_t skipped_peak_threshold = 0;
};

/**
 * @brief Apply snow-depth increments to Noah-MP snow state.
 *
 * This class is the C++ analogue of the Fortran NoahMPdisag_module.
 */
class NoahMPSnowIncrementer {
 public:
  /**
   * @param[in] vector_length Number of land-vector points.
   * @param[in] increment Snow-depth increments in millimeters.
   * @param[inout] state Noah-MP state over land points.
   * @param[in] noincr_threshold Skip positive increments when existing
   *                             snow depth exceeds this threshold.
   * @param[in] print_summary Print summary statistics.
   * @param[in] print_debug Reserved for compatibility with the Fortran API.
   * @return Summary statistics for the update.
   */
  static SnowIncrementSummary updateAllLayers(
      std::size_t vector_length,
      const std::vector<double>& increment,
      NoahMPSnowState& state,
      double noincr_threshold,
      bool print_summary,
      bool print_debug = false) {
    (void)print_debug;

    validateInputs(vector_length, increment, state);

    SnowIncrementSummary summary;

    for (std::size_t iloc = 0; iloc < vector_length; ++iloc) {
      const double temp_soil_corr =
          std::min(273.15, state.temperature_soil[iloc]);

      const double anal_snow_depth =
          state.snow_depth[iloc] + increment[iloc];

      // Same threshold as the Fortran implementation.
      if (std::abs(increment[iloc]) <= 0.01) {
        ++summary.no_increment;
        continue;
      }

      // The analysis has zero or negative snow depth.
      if (anal_snow_depth <= 0.0001) {
        clearSnowPoint(state, iloc);
        ++summary.removed_all_snow;
        continue;
      }

      const int active_layers =
          static_cast<int>(std::llround(state.active_snow_layers[iloc]));

      if (active_layers < 0) {
        updateMultilayerPoint(
            iloc,
            increment[iloc],
            active_layers,
            temp_soil_corr,
            noincr_threshold,
            state,
            summary);
      } else if (active_layers == 0) {
        updateZeroLayerPoint(
            iloc,
            increment[iloc],
            temp_soil_corr,
            state,
            summary);
      }
    }

    if (print_summary) {
      printSummary(summary);
    }

    return summary;
  }

 private:
  static void validateInputs(
      std::size_t vector_length,
      const std::vector<double>& increment,
      const NoahMPSnowState& state) {
    if (increment.size() < vector_length) {
      throw std::invalid_argument(
          "Snow increment vector is smaller than vector_length");
    }

    if (state.size() < vector_length) {
      throw std::invalid_argument(
          "Noah-MP state is smaller than vector_length");
    }

    if (!state.consistent()) {
      throw std::invalid_argument(
          "Noah-MP snow state arrays have inconsistent sizes");
    }
  }

  static void clearSnowPoint(
      NoahMPSnowState& state,
      std::size_t iloc) {
    state.swe[iloc] = 0.0;
    state.snow_depth[iloc] = 0.0;
    state.active_snow_layers[iloc] = 0.0;
    state.swe_previous[iloc] = 0.0;

    state.snow_soil_interface[iloc] = kZeroSnowInterfaces;
    state.temperature_snow[iloc].fill(0.0);
    state.snow_ice_layer[iloc].fill(0.0);
    state.snow_liq_layer[iloc].fill(0.0);
  }

  static void updateMultilayerPoint(
      std::size_t iloc,
      double increment,
      int active_layers,
      double /*temp_soil_corr*/,
      double noincr_threshold,
      NoahMPSnowState& state,
      SnowIncrementSummary& summary) {
    std::array<double, kSnowLayers> layer_depths{};

    layer_depths[0] =
        state.snow_soil_interface[iloc][0];

    layer_depths[1] =
        state.snow_soil_interface[iloc][1] -
        state.snow_soil_interface[iloc][0];

    layer_depths[2] =
        state.snow_soil_interface[iloc][2] -
        state.snow_soil_interface[iloc][1];

    if (increment > 0.0) {
      if (state.snow_depth[iloc] > noincr_threshold) {
        ++summary.skipped_peak_threshold;
      } else {
        ++summary.added_multilayer;

        // Fortran:
        //   vector_loc = 4 + active_layers
        //
        // The equivalent zero-based C++ layer index is:
        //   3 + active_layers
        const int first_layer = 3 + active_layers;

        updateLayers(
            iloc,
            first_layer,
            increment,
            layer_depths,
            state);
      }
    } else if (increment < 0.0) {
      ++summary.removed_multilayer;

      const int first_layer = 3 + active_layers;

      updateLayers(
          iloc,
          first_layer,
          increment,
          layer_depths,
          state);
    }

    // Recalculate interfaces 4..7 in Fortran, or 3..6 in C++.
    for (std::size_t ilayer = 3;
         ilayer < kSnowInterfaces;
         ++ilayer) {
      state.snow_soil_interface[iloc][ilayer] =
          state.snow_soil_interface[iloc][2] -
          kSoilInterfaces[ilayer];
    }

    state.snow_depth[iloc] =
        -state.snow_soil_interface[iloc][2] * 1000.0;

    state.swe[iloc] = 0.0;

    for (std::size_t ilayer = 0;
         ilayer < kSnowLayers;
         ++ilayer) {
      state.swe[iloc] +=
          state.snow_ice_layer[iloc][ilayer] +
          state.snow_liq_layer[iloc][ilayer];
    }

    state.swe_previous[iloc] = state.swe[iloc];

    // Exit multilayer mode below 25 mm.
    if (state.snow_depth[iloc] < 25.0) {
      state.active_snow_layers[iloc] = 0.0;
      state.snow_soil_interface[iloc] = kZeroSnowInterfaces;
      state.temperature_snow[iloc].fill(0.0);
      state.snow_ice_layer[iloc].fill(0.0);
      state.snow_liq_layer[iloc].fill(0.0);
    }
  }

  static void updateLayers(
      std::size_t iloc,
      int first_layer,
      double increment,
      const std::array<double, kSnowLayers>& layer_depths,
      NoahMPSnowState& state) {
    // The Fortran loop is from vector_loc through 3.
    // In C++, valid layer indices are 0 through 2.
    if (first_layer < 0) {
      first_layer = 0;
    }

    if (first_layer >= static_cast<int>(kSnowLayers)) {
      return;
    }

    for (int ilayer = first_layer;
         ilayer < static_cast<int>(kSnowLayers);
         ++ilayer) {
      const double depth = -layer_depths[ilayer];

      // This is normally positive for an active snow layer.
      if (depth <= std::numeric_limits<double>::epsilon()) {
        continue;
      }

      const double total_layer_water =
          state.snow_ice_layer[iloc][ilayer] +
          state.snow_liq_layer[iloc][ilayer];

      if (total_layer_water <= std::numeric_limits<double>::epsilon()) {
        continue;
      }

      const double partition_ratio =
          depth / state.snow_depth[iloc] * 1000.0;

      const double layer_density =
          total_layer_water / depth;

      const double swe_increment =
          partition_ratio * increment * layer_density / 1000.0;

      const double liq_ratio =
          state.snow_liq_layer[iloc][ilayer] /
          total_layer_water;

      state.snow_ice_layer[iloc][ilayer] +=
          (1.0 - liq_ratio) * swe_increment;

      state.snow_liq_layer[iloc][ilayer] +=
          liq_ratio * swe_increment;

      for (int iinter = ilayer;
           iinter < static_cast<int>(kSnowLayers);
           ++iinter) {
        state.snow_soil_interface[iloc][iinter] -=
            partition_ratio * increment / 1000.0;
      }
    }
  }

  static void updateZeroLayerPoint(
      std::size_t iloc,
      double increment,
      double temp_soil_corr,
      NoahMPSnowState& state,
      SnowIncrementSummary& summary) {
    if (increment > 0.0) {
      double layer_density = 0.0;

      if (state.snow_depth[iloc] < 1.0) {
        // New-snow density as used in the Fortran implementation.
        layer_density = std::max(
            80.0,
            std::min(
                120.0,
                67.92 +
                    51.25 *
                        std::exp((temp_soil_corr - 273.15) / 2.59)));
      } else {
        layer_density =
            state.swe[iloc] /
            state.snow_depth[iloc] *
            1000.0;
      }

      if (state.temperature_soil[iloc] <= 273.155) {
        const double delta =
            std::min(
                state.snow_depth[iloc] + increment,
                50.0) -
            state.snow_depth[iloc];

        state.snow_depth[iloc] =
            std::min(
                state.snow_depth[iloc] + increment,
                50.0);

        state.swe[iloc] +=
            delta * layer_density / 1000.0;

        ++summary.added_zero_layer;
      } else {
        ++summary.skipped_too_warm;
      }

      state.swe_previous[iloc] = state.swe[iloc];

      state.active_snow_layers[iloc] = 0.0;
      state.snow_ice_layer[iloc].fill(0.0);
      state.snow_liq_layer[iloc].fill(0.0);
      state.temperature_snow[iloc].fill(0.0);

      for (std::size_t ilayer = 0;
           ilayer < kSnowLayers;
           ++ilayer) {
        state.snow_soil_interface[iloc][ilayer] = 0.0;
      }

      // Promote zero-layer snow to multilayer mode.
      if (state.snow_depth[iloc] > 25.0) {
        ++summary.created_layer;

        state.active_snow_layers[iloc] = -1.0;
        state.snow_ice_layer[iloc][2] = state.swe[iloc];
        state.temperature_snow[iloc][2] = temp_soil_corr;

        for (std::size_t ilayer = 2;
             ilayer < kSnowInterfaces;
             ++ilayer) {
          state.snow_soil_interface[iloc][ilayer] -=
              state.snow_depth[iloc] / 1000.0;
        }
      }

    } else if (increment < 0.0) {
      ++summary.removed_zero_layer;

      // This follows the original Fortran implementation. The
      // zero-layer branch is expected to have positive snow depth.
      if (state.snow_depth[iloc] <= 0.0) {
        clearSnowPoint(state, iloc);
        return;
      }

      const double layer_density =
          state.swe[iloc] /
          state.snow_depth[iloc] *
          1000.0;

      state.snow_depth[iloc] += increment;
      state.swe[iloc] +=
          increment * layer_density / 1000.0;

      state.swe_previous[iloc] = state.swe[iloc];

      state.active_snow_layers[iloc] = 0.0;
      state.snow_ice_layer[iloc].fill(0.0);
      state.snow_liq_layer[iloc].fill(0.0);
      state.temperature_snow[iloc].fill(0.0);

      for (std::size_t ilayer = 0;
           ilayer < kSnowLayers;
           ++ilayer) {
        state.snow_soil_interface[iloc][ilayer] = 0.0;
      }
    }
  }

  static void printSummary(
      const SnowIncrementSummary& summary) {
    std::cout << "Noah-MP snow increment summary:\n";
    std::cout << "No increments added: "
              << summary.no_increment << '\n';
    std::cout << "Increment removed all snow: "
              << summary.removed_all_snow << '\n';
    std::cout << "Increment added snow in multi-layer mode: "
              << summary.added_multilayer << '\n';
    std::cout << "Increment removed snow in multi-layer mode: "
              << summary.removed_multilayer << '\n';
    std::cout << "Increment added snow in zero-layer mode: "
              << summary.added_zero_layer << '\n';
    std::cout << "Increment not added in zero-layer mode, too warm: "
              << summary.skipped_too_warm << '\n';
    std::cout << "Increment added in zero-layer mode, added a layer: "
              << summary.created_layer << '\n';
    std::cout << "Increment removed snow in zero-layer mode: "
              << summary.removed_zero_layer << '\n';
    std::cout << "Positive increment skipped on snow depth "
                 "exceeding peak threshold: "
              << summary.skipped_peak_threshold << '\n';
  }
};

}  // namespace snow_increments
