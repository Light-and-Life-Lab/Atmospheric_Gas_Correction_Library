#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

struct L1_Record_PY
{
public:
    L1_Record_PY(){};
    
    pybind11::array_t<double, pybind11::array::c_style> ozone_absorption_cross_section{};
    pybind11::array_t<double, pybind11::array::c_style> ozone_concentration{}; 
    pybind11::array_t<double, pybind11::array::c_style> no2_absorption_cross_section{};
    pybind11::array_t<double, pybind11::array::c_style> fraction_tropospheric_no2_above_200m{};
    pybind11::array_t<double, pybind11::array::c_style> tropospheric_no2_concentration{};
    pybind11::array_t<double, pybind11::array::c_style> stratospheric_no2_concentration{};
    pybind11::array_t<double, pybind11::array::c_style> co_transmittance{};
    pybind11::array_t<double, pybind11::array::c_style> air_mass_factor_mixed_gases{};
    int num_amf_grid_points{};
    pybind11::array_t<double, pybind11::array::c_style> cos_solar_zenith{}; 
    pybind11::array_t<double, pybind11::array::c_style> cos_sensor_zenith{};
};


struct Transmittance_Record_PY
{
public:
    Transmittance_Record_PY(){};
    
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_solar_zenith{};
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_sensor_zenith{};
    pybind11::array_t<double, pybind11::array::c_style> gas_transmittance_total{};
};