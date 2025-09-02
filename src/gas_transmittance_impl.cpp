#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <omp.h>

#include "gas_transmittance.h"

int32_t get_index_lowerbound(double* table_val, int num_val, float val) 
{
    int32_t index;
    int32_t i;

    for (i = 0; i < num_val; i++)
        if (val < table_val[i])
            break;
    index = std::max(i-1, 0);
    index = std::min(index, num_val-2);
    return index;
}

std::pair<int, double> get_amf_index_and_ratio(L1_Record* l1_rec, double amf_value)
{
    int index_amf = get_index_lowerbound(l1_rec->air_mass_factor_mixed_gases, l1_rec->num_amf_grid_points, amf_value);
    double ratio = (amf_value - l1_rec->air_mass_factor_mixed_gases[index_amf]) /
                    (l1_rec->air_mass_factor_mixed_gases[index_amf + 1] - l1_rec->air_mass_factor_mixed_gases[index_amf]);

    return std::make_pair(index_amf, ratio);
}


double interpolate_transmittance_to_amf(L1_Record* l1_rec, double* transmittance_table, int32_t index, double ratio)
{
    // In the case where amf correction is performed, the gas transmittance input file will 
    // have transmittance values as a 2D matrix. The transmittance values in this matrix are
    // a function of both wavelength and air mass factor (amf)

    // This desired amf value may fall between two of the points in the amf grid, so in general
    // it is necessary to interpolate the values in the transmittance table to obtain the 
    // transmittance at the desired amf value

    double transmittance_interpolated_to_amf = transmittance_table[index]*(1-ratio)
                                             + transmittance_table[index+1]*ratio;
    return transmittance_interpolated_to_amf;
}


int windex(float wave, double twave[], int ntwave) {
    int iw, index;
    double wdiff;
    double wdiffmin = 99999.;

    for (iw = 0; iw < ntwave; iw++) {

        /* break on exact match */
        if (twave[iw] == wave) {
            index = iw;
            break;
        }

        /* look for closest */
        wdiff = fabs(twave[iw] - wave);
        if (wdiff < wdiffmin) {
            wdiffmin = wdiff;
            index = iw;
        }
    }

    return (index);
}


float get_airmass_oxygen(L1_Record* l1rec, int32_t ip, double window1, double absorp_band, double window2)
{
    int32_t i;

    int32_t num_wavelengths = l1rec->num_wavelengths;
    int32_t row_offset = ip*num_wavelengths;
    double* wavelength_array = l1rec->wavelengths;
    double* Lt = &l1rec->Lt[row_offset];
    double u0 = l1rec->cos_solar_zenith[ip];
    double* F0 = l1rec->F0;
    double rhot[3];

    int band1 = windex(window1, wavelength_array, num_wavelengths);
    rhot[0] = M_PI*Lt[band1]/F0[band1]/u0;

    int band2 = windex(window2, wavelength_array, num_wavelengths);
    rhot[1] = M_PI*Lt[band2]/F0[band2]/u0;

    int band_absorp = windex(absorp_band, wavelength_array, num_wavelengths);
    rhot[2] = M_PI*Lt[band_absorp]/F0[band_absorp]/u0;

    double rhot_interp = rhot[0]+(absorp_band-window1)*(rhot[1]-rhot[0])/(window2-window1);

    // std::cout << "u0: " << u0 << '\n';

    // std::cout << "Lt[band1]: " << Lt[band1] << '\n';
    // std::cout << "Lt[band2]: " << Lt[band2] << '\n';
    // std::cout << "Lt[band_absorp]: " << Lt[band_absorp] << '\n';

    // std::cout << "F0[band1]: " << F0[band1] << '\n';
    // std::cout << "F0[band2]: " << F0[band2] << '\n';
    // std::cout << "F0[band_absorp]: " << F0[band_absorp] << '\n';

    // std::cout << "rhot[0]: " << rhot[0] << '\n';
    // std::cout << "rhot[1]: " << rhot[1] << '\n';
    // std::cout << "rhot[2]: " << rhot[2] << '\n';

    // std::cout << "rhot_interp: " << rhot_interp << '\n';

    double trans_o2_true = rhot[2]/rhot_interp;

    // std::cout << "trans_o2_true: " << trans_o2_true << '\n';

    int num_airmass = l1rec->num_amf_grid_points;
    int amf_table_row_offset = band_absorp*num_airmass;
    // std::cout << "amf_table_row_offset: " << amf_table_row_offset << '\n';
    for (i = 0; i < num_airmass; i++) 
    {
        if (trans_o2_true >= l1rec->o2_transmittance[amf_table_row_offset + i])
            break;
    }
    if (i == 0)
        i = 1;
    if (i == num_airmass)
        i = num_airmass - 1;

    // std::cout << "l1rec->air_mass_factor_mixed_gases[i]" << l1rec->air_mass_factor_mixed_gases[i] << '\n';
    // std::cout << "l1rec->air_mass_factor_mixed_gases[i] - l1rec->air_mass_factor_mixed_gases[i - 1]" << l1rec->air_mass_factor_mixed_gases[i] - l1rec->air_mass_factor_mixed_gases[i - 1] << '\n';

    // std::cout << "l1rec->o2_transmittance[amf_table_row_offset+i]" << l1rec->o2_transmittance[amf_table_row_offset+i] << '\n';
    // std::cout << "l1rec->o2_transmittance[amf_table_row_offset+i - 1]" << l1rec->o2_transmittance[amf_table_row_offset+i - 1] << '\n';
    // std::cout << "trans_o2_true - l1rec->o2_transmittance[amf_table_row_offset+i]" << trans_o2_true - l1rec->o2_transmittance[amf_table_row_offset+i] << '\n';
    // std::cout << "l1rec->o2_transmittance[amf_table_row_offset+i] - l1rec->o2_transmittance[amf_table_row_offset+i - 1]" << l1rec->o2_transmittance[amf_table_row_offset+i] - l1rec->o2_transmittance[amf_table_row_offset+i - 1] << '\n';

    double amf_interp = l1rec->air_mass_factor_mixed_gases[i] + (trans_o2_true - l1rec->o2_transmittance[amf_table_row_offset+i]) 
                 * (l1rec->air_mass_factor_mixed_gases[i] - l1rec->air_mass_factor_mixed_gases[i - 1]) 
                 / (l1rec->o2_transmittance[amf_table_row_offset+i] - l1rec->o2_transmittance[amf_table_row_offset+i - 1]);

    return (amf_interp);
}


void ozone_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            double ozone_optical_depth = l1_rec->ozone_concentration[ip] * l1_rec->ozone_absorption_cross_section[iw];
            t_rec->gas_transmittance_solar_zenith[row_offset + iw] = exp(-(ozone_optical_depth / l1_rec->cos_solar_zenith[ip]));

            if (do_amf_correction) 
            {
                t_rec->gas_transmittance_total[row_offset + iw] = exp(-ozone_optical_depth * (1.0/l1_rec->cos_solar_zenith[ip] + 1.0/l1_rec->cos_sensor_zenith[ip]));
            } 
            else 
            {
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = exp(-(ozone_optical_depth / l1_rec->cos_sensor_zenith[ip]));
            }
        }
    }
}


void co2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(l1_rec, amf_total);

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction)
            {
                int32_t row_index = iw*l1_rec->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co2_transmittance, table_index_solz, ratio_solz);
                t_rec->gas_transmittance_total[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co2_transmittance, table_index_total, ratio_total);
            }
            else
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->co2_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->co2_transmittance[iw], amf_sensor_zenith);
            }
        }
    }
}


void co_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction) 
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(l1_rec, amf_total);

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction)
            {
                int32_t row_index = iw*l1_rec->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co_transmittance, table_index_solz, ratio_solz);
                t_rec->gas_transmittance_total[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co_transmittance, table_index_total, ratio_total);
            }
            else
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->co_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->co_transmittance[iw], amf_sensor_zenith);
            }
        }
    }
}

void ch4_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(l1_rec, amf_total);

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction)
            {
                int32_t row_index = iw*l1_rec->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->ch4_transmittance, table_index_solz, ratio_solz);
                t_rec->gas_transmittance_total[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->ch4_transmittance, table_index_total, ratio_total);
            }
            else
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->ch4_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->ch4_transmittance[iw], amf_sensor_zenith);
            }
        }
    }
}

void o2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction, Oxygen_A_Band_Option oxygen_A_band_option) 
{
    // #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(l1_rec, amf_total);

        int index_amf_solz_o2;
        float ratio_solz_o2;
        float ratio_total_o2;
        int index_amf_total_o2;

        if (do_amf_correction && oxygen_A_band_option == Oxygen_A_Band_Option::YES_AMF_CORRECTION) 
        {
            float amf_total_o2 = get_airmass_oxygen(l1_rec, ip, 753.0221, 761.7891, 776.81335);
            float scaling_factor = amf_total_o2 / amf_total;

            std::cout << "amf_total_o2: " << amf_total_o2 << '\n';
            std::cout << "scaling_factor: " << scaling_factor << std::endl;

            // auto [index_amf_solz_o2, ratio_solz_o2] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith * scaling_factor);
            // auto [index_amf_total_o2, ratio_total_o2] = get_amf_index_and_ratio(l1_rec, amf_total * scaling_factor);

            index_amf_solz_o2 = get_index_lowerbound(l1_rec->air_mass_factor_mixed_gases, l1_rec->num_amf_grid_points, amf_solar_zenith * scaling_factor);
            index_amf_total_o2 = get_index_lowerbound(l1_rec->air_mass_factor_mixed_gases, l1_rec->num_amf_grid_points, amf_total * scaling_factor);

            ratio_solz_o2 = (amf_solar_zenith * scaling_factor - l1_rec->air_mass_factor_mixed_gases[index_amf_solz_o2]) /
                            (l1_rec->air_mass_factor_mixed_gases[index_amf_solz_o2 + 1] - l1_rec->air_mass_factor_mixed_gases[index_amf_solz_o2]);
            ratio_total_o2 = (amf_total * scaling_factor - l1_rec->air_mass_factor_mixed_gases[index_amf_total_o2]) /
                                (l1_rec->air_mass_factor_mixed_gases[index_amf_total_o2 + 1] - l1_rec->air_mass_factor_mixed_gases[index_amf_total_o2]);
        }

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction) 
            {
                int32_t index=iw*l1_rec->num_amf_grid_points;
                float t_o2_interp;

                if (oxygen_A_band_option == Oxygen_A_Band_Option::YES_AMF_CORRECTION) 
                {
                    t_o2_interp = l1_rec->o2_transmittance[index + index_amf_solz_o2] * (1 - ratio_solz_o2) +
                                    l1_rec->o2_transmittance[index + index_amf_solz_o2 + 1] * ratio_solz_o2;
                    t_rec->gas_transmittance_solar_zenith[row_offset + iw] = t_o2_interp;

                    t_o2_interp = l1_rec->o2_transmittance[index + index_amf_total_o2] * (1 - ratio_total_o2) +
                                    l1_rec->o2_transmittance[index + index_amf_total_o2 + 1] * ratio_total_o2;
                    t_rec->gas_transmittance_total[row_offset + iw] = t_o2_interp;
                }
                else
                {
                    t_o2_interp = l1_rec->o2_transmittance[index + index_amf_solz] * (1 - ratio_solz) +
                                    l1_rec->o2_transmittance[index + index_amf_solz + 1] * ratio_solz;
                    t_rec->gas_transmittance_solar_zenith[row_offset + iw] = t_o2_interp;

                    t_o2_interp = l1_rec->o2_transmittance[index + index_amf_total] * (1 - ratio_total) +
                                    l1_rec->o2_transmittance[index + index_amf_total + 1] * ratio_total;
                    t_rec->gas_transmittance_total[row_offset + iw] = t_o2_interp;
                }
            }
            else
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->o2_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->o2_transmittance[iw], amf_sensor_zenith);
            }
        }
    }
}


void n2o_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        double amf_total = amf_solar_zenith + amf_sensor_zenith;

        auto [index_amf_solz, ratio_solz] = get_amf_index_and_ratio(l1_rec, amf_solar_zenith);
        auto [index_amf_total, ratio_total] = get_amf_index_and_ratio(l1_rec, amf_total);

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction)
            {
                int32_t row_index = iw*l1_rec->num_amf_grid_points;
                int32_t table_index_solz = row_index + index_amf_solz;
                int32_t table_index_total = row_index + index_amf_total;

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->n2o_transmittance, table_index_solz, ratio_solz);
                t_rec->gas_transmittance_total[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->n2o_transmittance, table_index_total, ratio_total);
            }
            else
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->n2o_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->n2o_transmittance[iw], amf_sensor_zenith);
            }
        }
    }
}

void no2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        double sec0 = 1.0 / l1_rec->cos_solar_zenith[ip];
        double sec = 1.0 / l1_rec->cos_sensor_zenith[ip];
        double tropospheric_no2_concentration_above_200m{0.0};

        if (l1_rec->tropospheric_no2_concentration[ip] > 0.0)
        {
            /* compute tropo no2 above 200m (Z.Ahmad)
            tropospheric_no2_concentration_above_200m = exp(12.6615 + 0.61676*log(no2_tropo));
            new, location-dependent method */
            tropospheric_no2_concentration_above_200m = l1_rec->fraction_tropospheric_no2_above_200m[ip] * l1_rec->tropospheric_no2_concentration[ip];
        }

        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (l1_rec->no2_absorption_cross_section[iw] > 0.0) 
            {
                double a_285 = l1_rec->no2_absorption_cross_section[iw] * (1.0 - 0.003 * (285.0 - 294.0));
                double a_225 = l1_rec->no2_absorption_cross_section[iw] * (1.0 - 0.003 * (225.0 - 294.0));

                double no2_optical_depth_to_200m = a_285 * tropospheric_no2_concentration_above_200m 
                                                 + a_225 * l1_rec->stratospheric_no2_concentration[ip];

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = exp(-(no2_optical_depth_to_200m * sec0));

                if (do_amf_correction) 
                {
                    t_rec->gas_transmittance_total[row_offset + iw] = exp(-(no2_optical_depth_to_200m * (sec + sec0)));
                }
                else
                {
                    t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = exp(-(no2_optical_depth_to_200m * sec));
                }
            }
        }
    }
}