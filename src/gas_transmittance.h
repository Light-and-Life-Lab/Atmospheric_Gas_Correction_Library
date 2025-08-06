struct L1_Record
{
    double* k_oz{};
    double* l1b_oz{};

    double* k_no2{};
    double* l1b_no2_frac{};
    double* l1b_no2_tropo{};
    double* l1b_no2_strat{};

    double* t_co{};

    double* amf_mixed{};
    int num_airmass{};

    double* l1b_csolz{}; 
    double* l1b_csenz{};
    int num_pixels{};
    int num_wavelengths{};
};

struct Transmittance_Record
{
    double* tg_sol{};
    double* tg_sen{};
    double* tg{};
};

void ozone_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void no2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void co_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);