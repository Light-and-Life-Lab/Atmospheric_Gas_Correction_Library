#ifndef PYBIND_INTERFACE_TYPES_H
#define PYBIND_INTERFACE_TYPES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "gas_transmittance.h"

struct Ancillary_Data_PY
{
public:
    Ancillary_Data_PY(){};

    pybind11::array_t<double, pybind11::array::c_style> ozone_absorption_cross_section{};
    pybind11::array_t<double, pybind11::array::c_style> ozone_concentration{}; 

    pybind11::array_t<double, pybind11::array::c_style> no2_absorption_cross_section{};
    pybind11::array_t<double, pybind11::array::c_style> fraction_tropospheric_no2_above_200m{};
    pybind11::array_t<double, pybind11::array::c_style> tropospheric_no2_concentration{};
    pybind11::array_t<double, pybind11::array::c_style> stratospheric_no2_concentration{};

    pybind11::array_t<double, pybind11::array::c_style> precipitable_water{};

    pybind11::array_t<double, pybind11::array::c_style> a_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> b_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> c_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> d_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> e_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> f_h2o{};
    pybind11::array_t<double, pybind11::array::c_style> g_h2o{};

    pybind11::array_t<double, pybind11::array::c_style> water_vapor_bands{};
    int num_water_vapor_bands{};
};


struct Air_Mass_Factor_Lookup_Table_PY
{
public:
    Air_Mass_Factor_Lookup_Table_PY(){};
    pybind11::array_t<double, pybind11::array::c_style> co2_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> co_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> ch4_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> o2_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> n2o_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> h2o_transmittance{};

    int model{};
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_table_wavelengths{};
    pybind11::array_t<double, pybind11::array::c_style> air_mass_factor_mixed_gases{};
    pybind11::array_t<double, pybind11::array::c_style> air_mass_factor_water_vapor{};
    pybind11::array_t<double, pybind11::array::c_style> water_vapor_concentration{};
    int num_models{};
    int num_gas_transmittance_wavelengths{};
    int num_amf_grid_points{}; // Length of both mixed gases and water vapor air mass factor tables
    int num_water_vapor_concentrations{};
};

struct L1_Record_PY
{
public:
    L1_Record_PY(){};

    // TODO: Replace Lt and F0 with reflectance (rhot), but only once OCSSW is running and we can use unit tests to catch errors when refactoring
    pybind11::array_t<double, pybind11::array::c_style> Lt{};
    pybind11::array_t<double, pybind11::array::c_style> F0{};

    pybind11::array_t<double, pybind11::array::c_style> cos_solar_zenith{}; 
    pybind11::array_t<double, pybind11::array::c_style> cos_sensor_zenith{};

    int num_pixels{};
    int num_wavelengths{};

    pybind11::array_t<double, pybind11::array::c_style> wavelengths{};
};


struct Transmittance_Record_PY
{
public:
    Transmittance_Record_PY(){};
    
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_solar_zenith{};
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_sensor_zenith{};
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_total{};
};

#endif // PYBIND_INTERFACE_TYPES_H