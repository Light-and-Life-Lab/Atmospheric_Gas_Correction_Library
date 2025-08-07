struct L1_Record
{
    double* ozone_absorption_cross_section{};
    double* ozone_concentration{};

    double* no2_absorption_cross_section{};
    double* fraction_tropospheric_no2_above_200m{};
    double* tropospheric_no2_concentration{};
    double* stratospheric_no2_concentration{};

    double* co2_transmittance{};
    double* co_transmittance{};
    double* ch4_transmittance{};
    double* n2o_transmittance{};

    double* air_mass_factor_mixed_gases{};
    int num_amf_grid_points{};

    double* cos_solar_zenith{}; 
    double* cos_sensor_zenith{};
    int num_pixels{};
    int num_wavelengths{};
};

struct Transmittance_Record
{
    double* gas_transmittance_solar_zenith{};
    double* gas_transmittance_sensor_zenith{};
    double* gas_transmittance_total{};
};

void ozone_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void co2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void co_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void ch4_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void n2o_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);
void no2_transmittance(L1_Record* l1_rec, Transmittance_Record* t_rec, bool do_amf_correction);