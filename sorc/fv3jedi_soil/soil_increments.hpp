/**
 * @file soil_increments.hpp
 * @brief Routines for applying soil DA increments
 * Copied from land_increments
 * @author Clara Draper ESRL/PSL
 * @author Tseganeh ZG April 2026 - bring soil specific parts close to snow DA increments code in GDASApp
 */

#ifndef SOIL_INCREMENTS_HPP
#define SOIL_INCREMENTS_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>

namespace soil_increments {

// Constants
constexpr int LSM_NOAHMP = 2;           ///< Flag for NOAHMP land surface model
constexpr double TFREEZ = 273.16;       ///< con_t0c in physcons
constexpr double HFUS = 0.3336e06;      ///< Latent heat of fusion (J/kg)
constexpr double GRAV = 9.80616;        ///< Gravity acceleration (m/s^2)

// Constants for frh2o function
constexpr double GS2 = 9.81;            ///< con_g in snowpack, frh2o
constexpr double LSUBF = 3.335e5;       ///< con_hfus=3.3358e+5
constexpr double CK = 8.0;
constexpr double BLIM = 5.5;
constexpr double ERROR = 0.005;

/**
 * @brief Calculate the liquid water (slc) for a given total moisture content 
 * and soil temperature. Used here to update slc when DA update to stc crosses 
 * freezing.
 * 
 * Subroutine frh2o calculates amount of supercooled liquid soil water content 
 * if temperature is below 273.15K. Requires Newton-type iteration to solve 
 * the nonlinear implicit equation given in eqn 17 of Koren et al (1999, JGR, 
 * vol 104(d16), 19569-19585).
 *
 * @param[in] tkelv Soil temperature in K
 * @param[in] smc Total soil moisture content (volumetric)
 * @param[in] sh2o Liquid soil moisture content (volumetric)
 * @param[in] smcmax Saturation soil moisture content
 * @param[in] bexp Soil type "b" parameter
 * @param[in] psis Saturated soil matric potential
 * @param[out] liqwat Supercooled liquid water content
 *
 * @author From Noah LSM sflx.F
 */
void frh2o(
    double tkelv,
    double smc,
    double sh2o,
    double smcmax,
    double bexp,
    double psis,
    double& liqwat)
{
    double bx = bexp;
    if (bexp > BLIM) {
        bx = BLIM;
    }

    int nlog = 0;
    int kcount = 0;

    // If temperature not significantly below freezing, liqwat = smc
    if (tkelv > (TFREEZ - 1.e-3)) {
        liqwat = smc;
    } else {
        if (CK != 0.0) {
            // Option 1: Iterated solution for nonzero ck
            // Initial guess for swl (frozen content)
            double swl = smc - sh2o;

            // Keep within bounds
            swl = std::max(std::min(swl, smc - 0.02), 0.0);

            // Start iterations
            while ((nlog < 10) && (kcount == 0)) {
                nlog++;

                double df = std::log((psis * GS2 / LSUBF) * 
                           std::pow((1.0 + CK * swl), 2.0) *
                           std::pow((smcmax / (smc - swl)), bx)) - 
                           std::log(-(tkelv - TFREEZ) / tkelv);

                double denom = 2.0 * CK / (1.0 + CK * swl) + bx / (smc - swl);
                double swlk = swl - df / denom;

                // Bounds useful for mathematical solution
                swlk = std::max(std::min(swlk, smc - 0.02), 0.0);

                // Mathematical solution bounds applied
                double dswl = std::abs(swlk - swl);
                swl = swlk;

                // If dswl less or eq. error, no more iterations required
                if (dswl <= ERROR) {
                    kcount++;
                }
            }

            // Bounds applied are valid for physical solution
            liqwat = smc - swl;
        }

        // Option 2: Explicit solution for Flerchinger equation (ck=0)
        // Apply physical bounds to Flerchinger solution
        if (kcount == 0) {
            double fk = std::pow((LSUBF / (GS2 * (-psis))) * 
                       ((tkelv - TFREEZ) / tkelv), -1.0 / bx) * smcmax;

            fk = std::max(fk, 0.02);
            liqwat = std::min(fk, smc);
        }
    }
}

/**
 * @brief Initialize soil and vegetation parameters for Noah LSM.
 * This subroutine initializes soil and vegetation parameters needed 
 * for land increment calculations using Noah LSM.
 *
 * @param[in] isot Soil type index (currently only supports isot=1)
 * @param[in] ivet Vegetation type (not used for Noah LSM, kept for consistency)
 * @param[out] maxsmc Maximum soil moisture for each soil type (30 values)
 * @param[out] bb B exponent for each soil type (30 values)
 * @param[out] satpsi Saturated matric potential for each soil type (30 values)
 * @param[out] iret Return integer (0 on success, -1 on error)
 *
 * @author Clara Draper
 */
void set_soilveg_noah(
    int isot,
    int ivet,
    std::vector<double>& maxsmc,
    std::vector<double>& bb,
    std::vector<double>& satpsi,
    int& iret)
{
    maxsmc.resize(30);
    bb.resize(30);
    satpsi.resize(30);

    if (isot == 1) {
        // Using STASGO table
        bb = {4.05,  4.26, 4.74, 5.33, 5.33,  5.25,
              6.77,  8.72,  8.17, 10.73, 10.39,  11.55,
              5.25,  4.26,  4.05, 4.26,  11.55,  4.05,
              4.05,  0.00,  0.00, 0.00,  0.00,  0.00,
              0.00,  0.00,  0.00, 0.00,  0.00,  0.00};

        maxsmc = {0.395, 0.421, 0.434, 0.476, 0.476, 0.439,
                  0.404, 0.464, 0.465, 0.406, 0.468, 0.457,
                  0.464, 0.421, 0.200, 0.421, 0.457, 0.200,
                  0.395, 0.000, 0.000, 0.000, 0.000, 0.000,
                  0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

        satpsi = {0.035, 0.0363, 0.1413, 0.7586, 0.7586, 0.3548,
                  0.1349, 0.6166, 0.2630, 0.0977, 0.3236, 0.4677,
                  0.3548, 0.0363, 0.0350, 0.0363, 0.4677, 0.0350,
                  0.0350, 0.00, 0.00, 0.00, 0.00, 0.00,
                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

        iret = 0;
    } else {
        std::cerr << "set_soilveg_noah not coded for soil type " << isot << std::endl;
        iret = -1;
    }
}

/**
 * @brief Initialize soil and vegetation parameters for Noah-MP LSM.
 * This subroutine initializes soil and vegetation parameters needed 
 * for land increment calculations using Noah-MP LSM.
 * 
 * Noah-MP related parameters were extracted from noahmp_table.f.
 * isot (soil type) = 1: STATSGO must be selected if NoahMP is used.
 * ivet (vegetation type) = 1: IBGP is used by UFS offline Land DA for Noah-MP 
 * as of 07/13/2023.
 *
 * @param[in] isot Soil type index (currently only supports isot=1 for STATSGO)
 * @param[in] ivet Vegetation type (not used for now)
 * @param[out] maxsmc Maximum soil moisture for each soil type (30 values)
 * @param[out] bb B exponent for each soil type (30 values)
 * @param[out] satpsi Saturated matric potential for each soil type (30 values)
 * @param[out] iret Return integer (0 on success, -1 on error)
 *
 * @author Yuan Xue
 */
void set_soilveg_noahmp(
    int isot,
    int ivet,
    std::vector<double>& maxsmc,
    std::vector<double>& bb,
    std::vector<double>& satpsi,
    int& iret)
{
    maxsmc.resize(30);
    bb.resize(30);
    satpsi.resize(30);

    if (isot == 1) {
        // Set soil-dependent params (STATSGO is the only option for UFS, 07/13/2023)
        maxsmc = {0.339, 0.421, 0.434, 0.476, 0.484,
                  0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468,
                  0.439, 1.000, 0.200, 0.421, 0.468, 0.200,
                  0.339, 0.339, 0.000, 0.000, 0.000, 0.000,
                  0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

        bb = {2.79,  4.26, 4.74, 5.33, 3.86,  5.25,
              6.77,  8.72,  8.17, 10.73,  10.39, 11.55,
              5.25,  0.0,  2.79, 4.26,  11.55,  2.79,
              2.79,  0.00,  0.00, 0.00,  0.00,  0.00,
              0.00,  0.00,  0.00, 0.00,  0.00,  0.00};

        satpsi = {0.069, 0.036, 0.141, 0.759, 0.955,
                  0.355, 0.135, 0.617, 0.263, 0.098, 0.324, 0.468,
                  0.355, 0.00, 0.069, 0.036, 0.468, 0.069,
                  0.069, 0.00, 0.00, 0.00, 0.00, 0.00,
                  0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

        iret = 0;
    } else {
        std::cerr << "For Noah-MP, set_soilveg is not supported for soil type " << isot << std::endl;
        iret = -1;
    }
}

/**
 * @brief Read in soil state increments and add to soil states.
 * Adapted from original add_gsi_increment_soil routine.
 *
 * @param[in] lsoil_incr Number of soil layers (from top) to apply soil increments to
 * @param[in] stcinc Soil temperature increments on the cubed-sphere tile [lensfc, lsoil]
 * @param[in] slcinc Liquid soil moisture increments on the cubed-sphere tile [lensfc, lsoil]
 * @param[inout] stc_state Soil temperature state vector [lensfc, lsoil]
 * @param[inout] smc_state Soil moisture (liquid plus solid) state vector [lensfc, lsoil]
 * @param[inout] slc_state Liquid soil moisture state vector [lensfc, lsoil]
 * @param[out] stc_updated Integer to record whether STC in each grid cell was updated [lensfc]
 * @param[out] slc_updated Integer to record whether SLC in each grid cell was updated [lensfc]
 * @param[in] soilsnow_tile Land mask for increments on the cubed-sphere tile [lensfc]
 * @param[in] soilsnow_fg_tile First guess land mask for increments on the cubed-sphere tile [lensfc]
 * @param[in] lensfc Number of land points on a tile
 * @param[in] lsoil Number of soil layers
 * @param[in] myrank MPI rank number
 * @param[in] upd_stc Update soil temperature flag
 * @param[in] upd_slc Update soil moisture flag
 * @param[in] print_summary Print summary statistics
 * @param[in] print_debug Print debug statistics
 *
 * @author Yuan Xue. 11/2023
 * @author Tseganeh ZG. April 2024 (calling from jedi-apply_lnd_inc)
 */
void add_increment_soil(
    int lsoil_incr,
    const std::vector<std::vector<double>>& stcinc,
    const std::vector<std::vector<double>>& slcinc,
    std::vector<std::vector<double>>& stc_state,
    std::vector<std::vector<double>>& smc_state,
    std::vector<std::vector<double>>& slc_state,
    std::vector<int>& stc_updated,
    std::vector<int>& slc_updated,
    const std::vector<int>& soilsnow_tile,
    const std::vector<int>& soilsnow_fg_tile,
    int lensfc,
    int lsoil,
    int myrank,
    bool upd_stc,
    bool upd_slc,
    bool print_summary,
    bool print_debug)
{
    stc_updated.assign(lensfc, 0);
    slc_updated.assign(lensfc, 0);

    if (print_summary && myrank == 0) {
        std::cout << std::endl;
        std::cout << "add soil increments on cubed-sphere tiles" << std::endl;
        std::cout << "updating soil temps: " << upd_stc << std::endl;
        std::cout << "updating soil moisture: " << upd_slc << std::endl;
        std::cout << "adding to first " << lsoil_incr << " surface layers only" << std::endl;
    }

    // Initialize variables for count statistics
    int nother = 0;           // grid cells not land
    int nsnowupd = 0;         // grid cells with snow (temperature not yet updated)
    int nslcupd = 0;          // grid cells that are updated
    int nstcupd = 0;          // grid cells that are updated
    int nfrozen = 0;          // not updated as frozen soil
    int nfrozen_upd = 0;      // not updated as frozen soil

    // Main loop over grid cells
    for (int ij = 0; ij < lensfc; ++ij) {
        int mask_tile = soilsnow_tile[ij];
        int mask_fg_tile = soilsnow_fg_tile[ij];

        // mask: 1 - soil, 2 - snow, 0 - land-ice, -1 - not land
        if (mask_tile <= 0) {
            nother++;
            continue;
        }

        // Skip if snow is present before or after snow update
        if (mask_fg_tile == 2 || mask_tile == 2) {
            nsnowupd++;
            continue;
        }

        // Update soil temperature grid cells
        if (mask_tile == 1) {
            bool soil_freeze = false;
            bool soil_ice = false;

            for (int k = 0; k < lsoil_incr; ++k) {
                if (stc_state[ij][k] < TFREEZ) {
                    soil_freeze = true;
                }
                if (smc_state[ij][k] - slc_state[ij][k] > 0.001) {
                    soil_ice = true;
                }

                if (upd_stc) {
                    stc_state[ij][k] += stcinc[ij][k];  // TODO: do not add if < min_inc
                    if (k == 0) {
                        stc_updated[ij] = 1;
                        nstcupd++;
                    }
                }

                if (stc_state[ij][k] < TFREEZ && !soil_freeze && k == 0) {
                    nfrozen_upd++;
                }

                // Do not do updates if this layer or any above is frozen
                if (!soil_freeze && !soil_ice) {
                    if (upd_slc) {
                        if (k == 0) {
                            nslcupd++;
                            slc_updated[ij] = 1;
                        }
                        // Apply zero limit here (higher, model-specific limits are later)
                        slc_state[ij][k] = std::max(slc_state[ij][k] + slcinc[ij][k], 0.0);
                        smc_state[ij][k] = std::max(smc_state[ij][k] + slcinc[ij][k], 0.0);
                    }
                } else {
                    if (k == 0) {
                        nfrozen++;
                    }
                }
            }
        }
    }

    if (print_summary && myrank == 0) {
        std::cout << " statistics of grids number processed for rank : " << std::setw(2) << myrank << std::endl;
        std::cout << " soil grid total" << std::setw(8) << lensfc << std::endl;
        std::cout << " soil grid cells slc updated = " << std::setw(8) << nslcupd << std::endl;
        std::cout << " soil grid cells stc updated = " << std::setw(8) << nstcupd << std::endl;
        std::cout << " soil grid cells not updated, frozen = " << std::setw(8) << nfrozen << std::endl;
        std::cout << " soil grid cells update, became frozen = " << std::setw(8) << nfrozen_upd << std::endl;
        std::cout << " (not updated yet) snow grid cells = " << std::setw(8) << nsnowupd << std::endl;
        std::cout << " grid cells, without soil or snow = " << std::setw(8) << nother << std::endl;
    }

    if (print_debug) {
        std::cout << " statistics of grids number processed for rank : " << std::setw(2) << myrank << std::endl;
        std::cout << " soil grid total" << std::setw(8) << lensfc << std::endl;
        std::cout << " soil grid cells slc updated = " << std::setw(8) << nslcupd << std::endl;
        std::cout << " soil grid cells stc updated = " << std::setw(8) << nstcupd << std::endl;
        std::cout << " soil grid cells not updated, frozen = " << std::setw(8) << nfrozen << std::endl;
        std::cout << " soil grid cells update, became frozen = " << std::setw(8) << nfrozen_upd << std::endl;
        std::cout << " (not updated yet) snow grid cells = " << std::setw(8) << nsnowupd << std::endl;
        std::cout << " grid cells, without soil or snow = " << std::setw(8) << nother << std::endl;
    }
}

/**
 * @brief Calculate soil mask for land on model grid.
 * Output is 1 - soil, 2 - snow-covered, 0 - land ice, -1 - not land.
 *
 * @param[in] swe Model snow water equivalent [lensfc]
 * @param[in] vtype Model vegetation type [lensfc]
 * @param[in] stype Soil type [lensfc]
 * @param[in] lensfc Number of land points for this tile
 * @param[in] veg_type_landice Value of vegetation class that indicates land-ice
 * @param[out] mask Land mask for increments [lensfc]
 *
 * @author Clara Draper @date March 2021
 * @author Yuan Xue - introduce stype to make the mask calculation more generic
 */
void calculate_landinc_mask(
    const std::vector<double>& swe,
    const std::vector<int>& vtype,
    const std::vector<int>& stype,
    int lensfc,
    int veg_type_landice,
    std::vector<int>& mask)
{
    mask.assign(lensfc, -1);  // not land

    for (int i = 0; i < lensfc; ++i) {
        if (stype[i] > 0) {
            if (swe[i] > 0.001) {
                mask[i] = 2;  // snow covered land
            } else {
                mask[i] = 1;  // non-snow covered land
            }
        }
        if (vtype[i] == veg_type_landice) {
            mask[i] = 0;  // land-ice
        }
    }
}

/**
 * @brief Make adjustments to dependent variables after applying land increments.
 * These adjustments are model-dependent, and are currently only coded in full for Noah LSM.
 *
 * For Noah-MP, the adjustment scheme as of 11/09/2023:
 * Case 1: frozen => frozen, recalculate slc, smc remains
 * Case 2: unfrozen => frozen, recalculate slc, smc remains
 * Case 3: frozen => unfrozen, melt all soil ice (if any)
 * Case 4: unfrozen => unfrozen along with other cases, do nothing
 *
 * @param[in] lsoil_incr Number of soil layers (from top) to apply soil increments to
 * @param[in] isot Integer code for the soil type data set
 * @param[in] ivegsrc Integer code for the vegetation type data set
 * @param[in] lensfc Number of land points for this tile
 * @param[in] lsoil Number of soil layers
 * @param[in] isoiltype Soil types [lensfc]
 * @param[in] mask Mask indicating surface type [lensfc]
 * @param[in] stc_bck Background soil temperature states [lensfc, lsoil]
 * @param[inout] stc_adj Analysis soil temperature states [lensfc, lsoil]
 * @param[inout] smc_adj Analysis soil moisture states [lensfc, lsoil]
 * @param[inout] slc_adj Analysis liquid soil moisture states [lensfc, lsoil]
 * @param[in] stc_updated Integer to record whether STC in each grid cell was updated [lensfc]
 * @param[in] slc_updated Integer to record whether SLC in each grid cell was updated [lensfc]
 * @param[in] zsoil Depth of bottom of each soil layer [lsoil]
 * @param[in] upd_stc Update soil temperature flag
 * @param[in] upd_slc Update soil moisture flag
 * @param[in] myrank MPI rank number
 * @param[in] print_summary Print summary statistics
 * @param[in] print_debug Print debug statistics
 *
 * @author Clara Draper @date April 2021
 * @author Tseganeh ZG. April 2024 (calling from jedi-apply_lnd_inc)
 */
void apply_land_da_adjustments_soil(
    int lsoil_incr,
    int isot,
    int ivegsrc,
    int lensfc,
    int lsoil,
    const std::vector<int>& isoiltype,
    const std::vector<int>& mask,
    const std::vector<std::vector<double>>& stc_bck,
    std::vector<std::vector<double>>& stc_adj,
    std::vector<std::vector<double>>& smc_adj,
    std::vector<std::vector<double>>& slc_adj,
    const std::vector<int>& stc_updated,
    const std::vector<int>& slc_updated,
    const std::vector<double>& zsoil,
    bool upd_stc,
    bool upd_slc,
    int myrank,
    bool print_summary,
    bool print_debug)
{
    int n_stc = 0;
    int n_slc = 0;

    if (upd_stc) {
        std::vector<double> maxsmc;
        std::vector<double> bb;
        std::vector<double> satpsi;
        
        int iret = 0;
        set_soilveg_noahmp(isot, ivegsrc, maxsmc, bb, satpsi, iret);
        
        if (iret < 0) {
            std::cerr << "FATAL ERROR: problem in set_soilveg_noahmp" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 10);
        }

        for (int i = 0; i < lensfc; ++i) {
            if (stc_updated[i] == 1) {
                n_stc++;
                int soiltype = isoiltype[i];

                for (int l = 0; l < lsoil_incr; ++l) {
                    // Case 1 & 2: frozen or becoming frozen
                    if (stc_adj[i][l] < TFREEZ) {
                        // Recompute supercool liquid water, smc_adj remains unchanged
                        double smp = HFUS * (TFREEZ - stc_adj[i][l]) / (GRAV * stc_adj[i][l]);  // (m)
                        double slc_new = maxsmc[soiltype] * 
                                       std::pow(smp / satpsi[soiltype], -1.0 / bb[soiltype]);
                        slc_adj[i][l] = std::max(std::min(slc_new, smc_adj[i][l]), 0.0);
                    }
                    // Case 3: frozen => unfrozen, melt all soil ice
                    if (stc_adj[i][l] > TFREEZ) {  // do not rely on stc_bck
                        slc_adj[i][l] = smc_adj[i][l];
                    }
                }
            }
        }
    }

    if (upd_slc) {
        std::vector<double> dz(lsoil);
        dz[0] = -zsoil[0];
        for (int l = 1; l < lsoil; ++l) {
            dz[l] = -zsoil[l] + zsoil[l - 1];
        }

        std::cout << "Applying soil moisture mins " << std::endl;

        for (int i = 0; i < lensfc; ++i) {
            if (slc_updated[i] == 1) {
                n_slc++;
                // Apply SM bounds (later: add upper SMC limit)
                for (int l = 0; l < lsoil_incr; ++l) {
                    // Noah-MP minimum is 1 mm per layer (in SMC)
                    double min_val = 0.001 / dz[l];
                    slc_adj[i][l] = std::max(min_val, slc_adj[i][l]);
                    smc_adj[i][l] = std::max(min_val, smc_adj[i][l]);
                }
            }
        }
    }

    if (print_summary && myrank == 0) {
        std::cout << "statistics of grids number processed for rank : " << std::setw(2) << myrank << std::endl;
        std::cout << " soil grid total" << std::setw(8) << lensfc << std::endl;
        std::cout << " soil grid cells with slc update" << std::setw(8) << n_slc << std::endl;
        std::cout << " soil grid cells with stc update" << std::setw(8) << n_stc << std::endl;
    }

    if (print_debug) {
        std::cout << "statistics of grids number processed for rank : " << std::setw(2) << myrank << std::endl;
        std::cout << " soil grid total" << std::setw(8) << lensfc << std::endl;
        std::cout << " soil grid cells with slc update" << std::setw(8) << n_slc << std::endl;
        std::cout << " soil grid cells with stc update" << std::setw(8) << n_stc << std::endl;
    }
}

}  // namespace soil_increments

#endif  // SOIL_INCREMENTS_HPP
