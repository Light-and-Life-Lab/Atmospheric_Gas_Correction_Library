#ifndef PYBIND_INTERFACE_TYPES_H
#define PYBIND_INTERFACE_TYPES_H

#include <optional>
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

    pybind11::array get_solar_zenith() const
    {
        return solar_zenith_;
    }

    pybind11::array get_sensor_zenith() const
    {
        return sensor_zenith_;
    }

    pybind11::array get_total() const
    {
        // Since total is always derivable from solar and sensor zenith transmittances, no need to store the whole total array in memory, especially if caller doesn't end up using it.
        // Lazily evaluate the total array so that it is only computed and taking up memory if the caller specifically asks for it.
        // std::call_once guarantees the computation happens exactly once even if multiple threads call get_total() concurrently on the same instance.
+        std::call_once(*total_once_flag_, [this]()
+        {
+            // Can just import numpy from pybind here to do the array multiplication
+            static pybind11::object np_multiply = pybind11::module_::import("numpy").attr("multiply");  // Needs to be static here so that numpy isn't re-imported every time this branch runs.
+            total_ = np_multiply(solar_zenith_, sensor_zenith_);
+        });

        return *total_;
    }

    void set_solar_zenith(const pybind11::array& solar_zenith)
    {
        solar_zenith_ = solar_zenith;
    }

    void set_sensor_zenith(const pybind11::array& sensor_zenith)
    {
        sensor_zenith_ = sensor_zenith;
    }
    
private:
    pybind11::array solar_zenith_{};
    pybind11::array sensor_zenith_{};
    mutable std::optional<pybind11::array> total_{};
    mutable std::unique_ptr<std::once_flag> total_once_flag_{std::make_unique<std::once_flag>()};
};

#endif // PYBIND_INTERFACE_TYPES_H