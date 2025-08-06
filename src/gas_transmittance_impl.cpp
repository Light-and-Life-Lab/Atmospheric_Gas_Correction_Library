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

    int index_amf = get_index_lowerbound(l1_rec->amf_mixed, l1_rec->num_airmass, amf_value);
    double ratio = (amf_value - l1_rec->amf_mixed[index_amf]) /
                        (l1_rec->amf_mixed[index_amf + 1] - l1_rec->amf_mixed[index_amf]);

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
            double tau_oz = l1_rec->l1b_oz[ip] * l1_rec->k_oz[iw];
            t_rec->tg_sol[row_offset + iw] = exp(-(tau_oz / l1_rec->l1b_csolz[ip]));

            if (do_amf_correction) 
            {
                t_rec->tg[row_offset + iw] = exp(-tau_oz * (1.0/l1_rec->l1b_csolz[ip] + 1.0/l1_rec->l1b_csenz[ip]));
            } 
            else 
            {
                t_rec->tg_sen[row_offset + iw] = exp(-(tau_oz / l1_rec->l1b_csenz[ip]));
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

        double amf_solz = 1.0/l1_rec->l1b_csolz[ip];
        double amf_senz = 1.0/l1_rec->l1b_csenz[ip];
        
        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (do_amf_correction) 
            {
                double amf_total = amf_solz + amf_senz;
                int32_t index = iw*l1_rec->num_airmass;

                t_rec->tg_sol[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->t_co, index, amf_solz);
                t_rec->tg[row_offset + iw] = interpolate_transmittance_to_amf(l1_rec, l1_rec->t_co, index, amf_total);
            }
            else 
            {
                t_rec->tg_sol[row_offset + iw] = pow(l1_rec->t_co[iw], amf_solz);
                t_rec->tg_sen[row_offset + iw] = pow(l1_rec->t_co[iw], amf_senz);
            }
        }
    }
}


void no2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction)
{
    #pragma omp parallel for
    for (int ip = 0; ip < l1_rec->num_pixels; ip++)
    {
        double sec0 = 1.0 / l1_rec->l1b_csolz[ip];
        double sec = 1.0 / l1_rec->l1b_csenz[ip];
        double no2_tr200{0.0};

        if (l1_rec->l1b_no2_tropo[ip] > 0.0)
        {
            /* compute tropo no2 above 200m (Z.Ahmad)
            no2_tr200 = exp(12.6615 + 0.61676*log(no2_tropo));
            new, location-dependent method */
            no2_tr200 = l1_rec->l1b_no2_frac[ip] * l1_rec->l1b_no2_tropo[ip];
        }

        int row_offset = ip*l1_rec->num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < l1_rec->num_wavelengths; iw++) 
        {
            if (l1_rec->k_no2[iw] > 0.0) 
            {
                double a_285 = l1_rec->k_no2[iw] * (1.0 - 0.003 * (285.0 - 294.0));
                double a_225 = l1_rec->k_no2[iw] * (1.0 - 0.003 * (225.0 - 294.0));

                double tau_to200 = a_285 * no2_tr200 + a_225 * l1_rec->l1b_no2_strat[ip];

                t_rec->tg_sol[row_offset + iw] = exp(-(tau_to200 * sec0));

                if (do_amf_correction) 
                {
                    t_rec->tg[row_offset + iw] = exp(-(tau_to200 * (sec + sec0)));
                }
                else
                {
                    t_rec->tg_sen[row_offset + iw] = exp(-(tau_to200 * sec));
                }
            }
        }
    }
}