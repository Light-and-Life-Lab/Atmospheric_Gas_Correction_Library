#ifndef PYBIND_INTERFACE_TYPES_H
#define PYBIND_INTERFACE_TYPES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "gas_transmittance.h"

struct Ancillary_Data_PY
{
public:
    Ancillary_Data_PY(){};

    pybind11::array ozone_absorption_cross_section{};
    pybind11::array ozone_concentration{}; 

    pybind11::array no2_absorption_cross_section{};
    pybind11::array fraction_tropospheric_no2_above_200m{};
    pybind11::array tropospheric_no2_concentration{};
    pybind11::array stratospheric_no2_concentration{};

    pybind11::array precipitable_water{};

    pybind11::array water_vapor_bands{};
    int num_water_vapor_bands{};
};


struct Gas_Transmittance_Lookup_Table_PY
{
public:
    Gas_Transmittance_Lookup_Table_PY(){};
    pybind11::array_t<double, pybind11::array::c_style> co2_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> co_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> ch4_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> o2_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> n2o_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> h2o_transmittance{};

    int model{};
    pybind11::array_t<double, pybind11::array::c_style> wavelengths{};
    pybind11::array_t<double, pybind11::array::c_style> air_mass_factor_mixed_gases{};
    pybind11::array_t<double, pybind11::array::c_style> air_mass_factor_water_vapor{};
    pybind11::array_t<double, pybind11::array::c_style> water_vapor_concentration{};
    int num_models{};
    int num_wavelengths{};
    int num_amf_grid_points{}; // Length of both mixed gases and water vapor air mass factor tables
    int num_water_vapor_concentrations{};
};

struct L1_Data_PY
{
public:
    L1_Data_PY(){};

    pybind11::array reflectance{};

    pybind11::array cos_solar_zenith{}; 
    pybind11::array cos_sensor_zenith{};

    pybind11::array latitude{};
    pybind11::array longitude{};

    pybind11::array wavelengths{};

    int num_pixels{};
    int num_wavelengths{};

};


struct Gas_Transmittances_PY
{
public:
    Gas_Transmittances_PY(){};
    
    pybind11::array solar_zenith{};
    pybind11::array sensor_zenith{};
    pybind11::array total{};
};

#endif // PYBIND_INTERFACE_TYPES_H