#include <cmath>
#include <algorithm>
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


double interpolate_transmittance_to_amf(L1_Record* l1_rec, double* transmittance_table, int32_t index, double amf_value)
{
    // In the case where amf correction is performed, the gas transmittance input file will 
    // have transmittance values as a 2D matrix. The transmittance values in this matrix are
    // a function of both wavelength and air mass factor (amf)

    // This desired amf value may fall between two of the points in the amf grid, so in general
    // it is necessary to interpolate the values in the transmittance table to obtain the 
    // transmittance at the desired amf value

    int index_amf = get_index_lowerbound(l1_rec->air_mass_factor_mixed_gases, l1_rec->num_amf_grid_points, amf_value);
    double ratio = (amf_value - l1_rec->air_mass_factor_mixed_gases[index_amf]) /
                        (l1_rec->air_mass_factor_mixed_gases[index_amf + 1] - l1_rec->air_mass_factor_mixed_gases[index_amf]);

    double transmittance_interpolated_to_amf = transmittance_table[index+index_amf]*(1-ratio)
                                             + transmittance_table[index+index_amf+1]*ratio;
    return transmittance_interpolated_to_amf;
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


void co_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction) 
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        double amf_solar_zenith = 1.0/l1_rec->cos_solar_zenith[ip];
        double amf_sensor_zenith = 1.0/l1_rec->cos_sensor_zenith[ip];
        
        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction) 
            {
                double amf_total = amf_solar_zenith + amf_sensor_zenith;
                int32_t row_index = iw*l1_rec->num_amf_grid_points;

                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co_transmittance, row_index, amf_solar_zenith);
                t_rec->gas_transmittance_total[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->co_transmittance, row_index, amf_total);
            }
            else 
            {
                t_rec->gas_transmittance_solar_zenith[row_offset + iw] = pow(l1_rec->co_transmittance[iw], amf_solar_zenith);
                t_rec->gas_transmittance_sensor_zenith[row_offset + iw] = pow(l1_rec->co_transmittance[iw], amf_sensor_zenith);
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