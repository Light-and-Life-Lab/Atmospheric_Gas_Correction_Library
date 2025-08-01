#include <cmath>
#include <omp.h>

#include "gas_transmittance.h"

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