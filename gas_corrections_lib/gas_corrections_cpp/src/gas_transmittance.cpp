#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "allocate_output_array.hpp"
#include "pybind_interface_types.hpp"
#include "gas_transmittance.h"

namespace py = pybind11;

static bool is_float32(const py::array& arr)
{
    // Returns true as long as the elements of the array are float32
    return arr.dtype().is(py::dtype::of<float>());
}

struct Named_Array
{
    const char* name;
    const py::array* array_pointer;   // Specifically use a pointer here rather than reference so that it is trivially copyable inside an initializer_list
};

static void warn_non_float32(const char* function_name, const char* field_name, const py::array& arr)
{
    std::string dtype_str = py::str(arr.dtype()).cast<std::string>();
    std::string message = "\n\033[33mFloat32Warning (" + std::string(function_name) + "): '" + field_name + "' is dtype " + dtype_str + " (not float32).\n\033[0m";
    PyErr_WarnEx(PyExc_RuntimeWarning, message.c_str(), 1);
}

static bool all_float32(const char* function_name, std::initializer_list<Named_Array> arrays)
{
    // Check each array's floating point type. Only returns true if all passed arrays have type float32. 
    // Issues a warning if any are found with float64 to let the user know that intended memory savings from using float32 will not apply unless ALL arrays 
    bool result = true;
    for (const auto& named_array : arrays)
    {
        if (!is_float32(*named_array.array_pointer))
        {
            warn_non_float32(function_name, named_array.name, *named_array.array_pointer);
            result = false;
        }
    }

    if (result == false)
    {
        std::string message = "\n\033[33mWarning: At least one float64 array was passed to " + std::string(function_name) + ". This function will fall back to float64"
        + " to prevent unintended errors due to loss of precision. If you were trying to save memory by using float32 arrays, resolve all Float32Warnings. "
        + "Only once all errors are resolved will " + std::string(function_name) + " use float32 arrays for caclulation.\n\033[0m";
        PyErr_WarnEx(PyExc_RuntimeWarning, message.c_str(), 1);
    }

    return result;
}

template<typename T>
static py::array_t<T, py::array::c_style | py::array::forcecast> pybind11_c_style_cast(const py::array& arr)
{
    // Arrays are accepted through the pybind interface using the dynamic pybind11::array type. 
    // We need to enforce C style memory layout, so dynamic cast to a C style array of the appropriate datatype T (either float or double).
    return arr.cast<py::array_t<T, py::array::c_style | py::array::forcecast>>();
}


template<typename T>
Gas_Transmittances_PY ozone_transmittance_py(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data) 
{
    auto ozone_absorption_cross_section_t = pybind11_c_style_cast<T>(ancillary_data.ozone_absorption_cross_section);
    auto ozone_concentration_t = pybind11_c_style_cast<T>(ancillary_data.ozone_concentration);

    Ancillary_Data<T> ancillary_data_c{};
    ancillary_data_c.ozone_absorption_cross_section = static_cast<T*>(ozone_absorption_cross_section_t.request().ptr);
    ancillary_data_c.ozone_concentration = static_cast<T*>(ozone_concentration_t.request().ptr);
    
    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);

    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    ozone_transmittance<T>(&l1_data_c, &ancillary_data_c, &gas_transmittances_c);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY ozone_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data)
{
    // Note: ancillary_data_c.ozone_absorption_cross_section is intentionally excluded from checking because 
    // the cross sections are only defined for at most a few hundred wavelengths, so the size of that array is not impactful enough to warrant worrying about float32 vs float64
    bool all_arrays_are_float32 = all_float32(
        "ozone_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
            {"ancillary_data.ozone_concentration", &ancillary_data.ozone_concentration}
        }
    );

    if (all_arrays_are_float32)
        return ozone_transmittance_py<float>(l1_data, ancillary_data);
    else
        return ozone_transmittance_py<double>(l1_data, ancillary_data);
}


template<typename T>
Gas_Transmittances_PY co2_transmittance_py(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.co2_transmittance = static_cast<double*>(gas_transmittance_table.co2_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);    

    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    co2_transmittance<T>(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY co2_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    bool all_arrays_are_float32 = all_float32(
        "co2_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
        }
    );

    if (all_arrays_are_float32)
        return co2_transmittance_py<float>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
    else
        return co2_transmittance_py<double>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
}


template<typename T>
Gas_Transmittances_PY co_transmittance_py(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.co_transmittance = static_cast<double*>(gas_transmittance_table.co_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);    
    
    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    co_transmittance<T>(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY co_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    bool all_arrays_are_float32 = all_float32(
        "co_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
        }
    );

    if (all_arrays_are_float32)
        return co_transmittance_py<float>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
    else
        return co_transmittance_py<double>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
}


template<typename T>
Gas_Transmittances_PY ch4_transmittance_py(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.ch4_transmittance = static_cast<double*>(gas_transmittance_table.ch4_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);    

    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    ch4_transmittance<T>(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY ch4_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    bool all_arrays_are_float32 = all_float32(
        "ch4_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
        }
    );

    if (all_arrays_are_float32)
        return ch4_transmittance_py<float>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
    else
        return ch4_transmittance_py<double>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
}


template<typename T>
Gas_Transmittances_PY o2_transmittance_py(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension, Oxygen_A_Band_Option oxygen_A_band_option) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.o2_transmittance = static_cast<double*>(gas_transmittance_table.o2_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    
    auto reflectance_t = pybind11_c_style_cast<T>(l1_data.reflectance);
    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);    
    auto wavelengths_t = pybind11_c_style_cast<T>(l1_data.wavelengths);

    L1_Data<T> l1_data_c{};
    l1_data_c.reflectance = static_cast<T*>(reflectance_t.request().ptr);
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;
    l1_data_c.wavelengths = static_cast<T*>(wavelengths_t.request().ptr);

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    o2_transmittance<T>(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension, oxygen_A_band_option);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY o2_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension, Oxygen_A_Band_Option oxygen_A_band_option) 
{
    // Note: l1_data_c.wavelengths is intentionally excluded from checking because there will be at most a few hundred wavelengths, 
    // so the size of that array is not impactful enough to warrant worrying about float32 vs float64.
    bool all_arrays_are_float32 = all_float32(
        "o2_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
            {"l1_data.reflectance", &l1_data.reflectance},
        }
    );

    if (all_arrays_are_float32)
        return o2_transmittance_py<float>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension, oxygen_A_band_option);
    else
        return o2_transmittance_py<double>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension, oxygen_A_band_option);
}


template<typename T>
Gas_Transmittances_PY n2o_transmittance_py(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.n2o_transmittance = static_cast<double*>(gas_transmittance_table.n2o_transmittance.request().ptr);
    gas_transmittance_table_c.air_mass_factor_mixed_gases = static_cast<double*>(gas_transmittance_table.air_mass_factor_mixed_gases.request().ptr);
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);    
    
    L1_Data<T> l1_data_c{};

    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    n2o_transmittance<T>(&l1_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY n2o_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    bool all_arrays_are_float32 = all_float32(
        "n2o_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
        }
    );

    if (all_arrays_are_float32)
        return n2o_transmittance_py<float>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
    else
        return n2o_transmittance_py<double>(l1_data, gas_transmittance_table, lookup_table_has_amf_dimension);
}


template<typename T>
Gas_Transmittances_PY no2_transmittance_py(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data)
{
    auto no2_absorption_cross_section_t = pybind11_c_style_cast<T>(ancillary_data.no2_absorption_cross_section);
    auto fraction_tropospheric_no2_above_200m_t = pybind11_c_style_cast<T>(ancillary_data.fraction_tropospheric_no2_above_200m);
    auto tropospheric_no2_concentration_t = pybind11_c_style_cast<T>(ancillary_data.tropospheric_no2_concentration);
    auto stratospheric_no2_concentration_t = pybind11_c_style_cast<T>(ancillary_data.stratospheric_no2_concentration);

    Ancillary_Data<T> ancillary_data_c{};
    ancillary_data_c.no2_absorption_cross_section = static_cast<T*>(no2_absorption_cross_section_t.request().ptr);
    ancillary_data_c.fraction_tropospheric_no2_above_200m = static_cast<T*>(fraction_tropospheric_no2_above_200m_t.request().ptr);
    ancillary_data_c.tropospheric_no2_concentration = static_cast<T*>(tropospheric_no2_concentration_t.request().ptr);
    ancillary_data_c.stratospheric_no2_concentration = static_cast<T*>(stratospheric_no2_concentration_t.request().ptr);

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith); 
    
    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};

    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    no2_transmittance<T>(&l1_data_c, &ancillary_data_c, &gas_transmittances_c);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY no2_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data)
{
    // Note: ancillary_data_c.no2_absorption_cross_section is intentionally excluded from checking because 
    // the cross sections are only defined for at most a few hundred wavelengths, so the size of that array is not impactful enough to warrant worrying about float32 vs float64
    bool all_arrays_are_float32 = all_float32(
        "no2_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
            {"ancillary_data.fraction_tropospheric_no2_above_200m", &ancillary_data.fraction_tropospheric_no2_above_200m},
            {"ancillary_data.tropospheric_no2_concentration", &ancillary_data.tropospheric_no2_concentration},
            {"ancillary_data.stratospheric_no2_concentration", &ancillary_data.stratospheric_no2_concentration},
        }
    );

    if (all_arrays_are_float32)
        return no2_transmittance_py<float>(l1_data, ancillary_data);
    else
        return no2_transmittance_py<double>(l1_data, ancillary_data);
}


template<typename T>
Gas_Transmittances_PY h2o_transmittance_py(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    auto precipitable_water_t = pybind11_c_style_cast<T>(ancillary_data.precipitable_water);
    auto water_vapor_bands_t = pybind11_c_style_cast<T>(ancillary_data.water_vapor_bands);

    Ancillary_Data<T> ancillary_data_c{};
    ancillary_data_c.precipitable_water = static_cast<T*>(precipitable_water_t.request().ptr);
    ancillary_data_c.water_vapor_bands = static_cast<T*>(water_vapor_bands_t.request().ptr);
    ancillary_data_c.num_water_vapor_bands = ancillary_data.num_water_vapor_bands;

    Gas_Transmittance_Lookup_Table gas_transmittance_table_c{};
    gas_transmittance_table_c.h2o_transmittance = static_cast<double*>(gas_transmittance_table.h2o_transmittance.request().ptr);
    gas_transmittance_table_c.model = gas_transmittance_table.model;
    gas_transmittance_table_c.air_mass_factor_water_vapor = static_cast<double*>(gas_transmittance_table.air_mass_factor_water_vapor.request().ptr);
    gas_transmittance_table_c.wavelengths = static_cast<double*>(gas_transmittance_table.wavelengths.request().ptr);
    gas_transmittance_table_c.water_vapor_concentration = static_cast<double*>(gas_transmittance_table.water_vapor_concentration.request().ptr);
    gas_transmittance_table_c.num_models = gas_transmittance_table.num_models;
    gas_transmittance_table_c.num_wavelengths = gas_transmittance_table.num_wavelengths;
    gas_transmittance_table_c.num_amf_grid_points = gas_transmittance_table.num_amf_grid_points;
    gas_transmittance_table_c.num_water_vapor_concentrations = gas_transmittance_table.num_water_vapor_concentrations;

    auto cos_solar_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_solar_zenith);
    auto cos_sensor_zenith_t = pybind11_c_style_cast<T>(l1_data.cos_sensor_zenith);
    auto reflectance_t = pybind11_c_style_cast<T>(l1_data.reflectance);
    auto wavelengths_t = pybind11_c_style_cast<T>(l1_data.wavelengths);   
    
    L1_Data<T> l1_data_c{};
    l1_data_c.cos_solar_zenith = static_cast<T*>(cos_solar_zenith_t.request().ptr);
    l1_data_c.cos_sensor_zenith = static_cast<T*>(cos_sensor_zenith_t.request().ptr);
    l1_data_c.reflectance = static_cast<T*>(reflectance_t.request().ptr);
    l1_data_c.wavelengths = static_cast<T*>(wavelengths_t.request().ptr);
    l1_data_c.num_pixels = l1_data.num_pixels;
    l1_data_c.num_wavelengths = l1_data.num_wavelengths;

    int n_rows = l1_data_c.num_pixels;
    int n_cols = l1_data_c.num_wavelengths;

    auto solar_zenith = allocate_output_array<T>(n_rows, n_cols);
    auto sensor_zenith = allocate_output_array<T>(n_rows, n_cols);

    Gas_Transmittances<T> gas_transmittances_c{};
    gas_transmittances_c.solar_zenith = static_cast<T*>(solar_zenith.request().ptr);
    gas_transmittances_c.sensor_zenith = static_cast<T*>(sensor_zenith.request().ptr);

    h2o_transmittance<T>(&l1_data_c, &ancillary_data_c, &gas_transmittance_table_c, &gas_transmittances_c, lookup_table_has_amf_dimension);

    Gas_Transmittances_PY gas_transmittances{};
    gas_transmittances.set_solar_zenith(solar_zenith);
    gas_transmittances.set_sensor_zenith(sensor_zenith);

    return gas_transmittances;
}

Gas_Transmittances_PY h2o_transmittance_type_dispatcher(const L1_Data_PY& l1_data, const Ancillary_Data_PY& ancillary_data, const Gas_Transmittance_Lookup_Table_PY& gas_transmittance_table, const bool lookup_table_has_amf_dimension) 
{
    // Note: l1_data_c.wavelengths is intentionally excluded from checking because there will be at most a few hundred wavelengths, 
    // so the size of that array is not impactful enough to warrant worrying about float32 vs float64.
    // Similarly, ancillary_data_c.water_vapor_bands is only 3 values, not worth worrying about.
    bool all_arrays_are_float32 = all_float32(
        "h2o_transmittance", 
        {
            {"l1_data.cos_solar_zenith", &l1_data.cos_solar_zenith},
            {"l1_data.cos_sensor_zenith", &l1_data.cos_sensor_zenith},
            {"l1_data.reflectance", &l1_data.reflectance},
            {"ancillary_data.precipitable_water", &ancillary_data.precipitable_water},
        }
    );

    if (all_arrays_are_float32)
        return h2o_transmittance_py<float>(l1_data, ancillary_data, gas_transmittance_table, lookup_table_has_amf_dimension);
    else
        return h2o_transmittance_py<double>(l1_data, ancillary_data, gas_transmittance_table, lookup_table_has_amf_dimension);
}


PYBIND11_MODULE(gas_transmittance, m) 
{
    py::class_<Ancillary_Data_PY>(m, "Ancillary_Data", py::module_local())
        .def(py::init<>())
        .def_readwrite("ozone_absorption_cross_section", &Ancillary_Data_PY::ozone_absorption_cross_section)
        .def_readwrite("ozone_concentration", &Ancillary_Data_PY::ozone_concentration)
        .def_readwrite("no2_absorption_cross_section", &Ancillary_Data_PY::no2_absorption_cross_section)
        .def_readwrite("fraction_tropospheric_no2_above_200m", &Ancillary_Data_PY::fraction_tropospheric_no2_above_200m)
        .def_readwrite("tropospheric_no2_concentration", &Ancillary_Data_PY::tropospheric_no2_concentration)
        .def_readwrite("stratospheric_no2_concentration", &Ancillary_Data_PY::stratospheric_no2_concentration)
        .def_readwrite("precipitable_water", &Ancillary_Data_PY::precipitable_water)
        .def_readwrite("water_vapor_bands", &Ancillary_Data_PY::water_vapor_bands)
        .def_readwrite("num_water_vapor_bands", &Ancillary_Data_PY::num_water_vapor_bands);

    py::class_<Gas_Transmittance_Lookup_Table_PY>(m, "Gas_Transmittance_Lookup_Table", py::module_local())
        .def(py::init<>())
        .def_readwrite("co2_transmittance", &Gas_Transmittance_Lookup_Table_PY::co2_transmittance)
        .def_readwrite("co_transmittance", &Gas_Transmittance_Lookup_Table_PY::co_transmittance)
        .def_readwrite("ch4_transmittance", &Gas_Transmittance_Lookup_Table_PY::ch4_transmittance)
        .def_readwrite("o2_transmittance", &Gas_Transmittance_Lookup_Table_PY::o2_transmittance)
        .def_readwrite("n2o_transmittance", &Gas_Transmittance_Lookup_Table_PY::n2o_transmittance)
        .def_readwrite("h2o_transmittance", &Gas_Transmittance_Lookup_Table_PY::h2o_transmittance)
        .def_readwrite("model", &Gas_Transmittance_Lookup_Table_PY::model)
        .def_readwrite("wavelengths", &Gas_Transmittance_Lookup_Table_PY::wavelengths)
        .def_readwrite("air_mass_factor_mixed_gases", &Gas_Transmittance_Lookup_Table_PY::air_mass_factor_mixed_gases)
        .def_readwrite("air_mass_factor_water_vapor", &Gas_Transmittance_Lookup_Table_PY::air_mass_factor_water_vapor)
        .def_readwrite("water_vapor_concentration", &Gas_Transmittance_Lookup_Table_PY::water_vapor_concentration)
        .def_readwrite("num_models", &Gas_Transmittance_Lookup_Table_PY::num_models)
        .def_readwrite("num_wavelengths", &Gas_Transmittance_Lookup_Table_PY::num_wavelengths)
        .def_readwrite("num_amf_grid_points", &Gas_Transmittance_Lookup_Table_PY::num_amf_grid_points)
        .def_readwrite("num_water_vapor_concentrations", &Gas_Transmittance_Lookup_Table_PY::num_water_vapor_concentrations);

    py::class_<L1_Data_PY>(m, "L1_Data", py::module_local())
        .def(py::init<>())
        .def_readwrite("reflectance", &L1_Data_PY::reflectance)
        .def_readwrite("cos_solar_zenith", &L1_Data_PY::cos_solar_zenith)
        .def_readwrite("cos_sensor_zenith", &L1_Data_PY::cos_sensor_zenith)
        .def_readwrite("latitude", &L1_Data_PY::latitude)
        .def_readwrite("longitude", &L1_Data_PY::longitude)
        .def_readwrite("num_pixels", &L1_Data_PY::num_pixels)
        .def_readwrite("num_wavelengths", &L1_Data_PY::num_wavelengths)
        .def_readwrite("wavelengths", &L1_Data_PY::wavelengths);

    py::class_<Gas_Transmittances_PY>(m, "Gas_Transmittances", py::module_local())
        .def(py::init<>())
        .def_property_readonly("solar_zenith", &Gas_Transmittances_PY::get_solar_zenith)
        .def_property_readonly("sensor_zenith", &Gas_Transmittances_PY::get_sensor_zenith)
        .def_property_readonly("total", &Gas_Transmittances_PY::get_total);

    py::enum_<Oxygen_A_Band_Option>(m, "Oxygen_A_Band_Option", py::module_local())
        .value("TRANSMITTANCE_TABLE", Oxygen_A_Band_Option::TRANSMITTANCE_TABLE)
        .value("SURROUNDING_WINDOW_BANDS", Oxygen_A_Band_Option::SURROUNDING_WINDOW_BANDS);

    m.def("ozone_transmittance", &ozone_transmittance_type_dispatcher);
    m.def("co2_transmittance", &co2_transmittance_type_dispatcher);
    m.def("co_transmittance", &co_transmittance_type_dispatcher);
    m.def("ch4_transmittance", &ch4_transmittance_type_dispatcher);
    m.def("o2_transmittance", &o2_transmittance_type_dispatcher);
    m.def("n2o_transmittance", &n2o_transmittance_type_dispatcher);
    m.def("no2_transmittance", &no2_transmittance_type_dispatcher);
    m.def("h2o_transmittance", &h2o_transmittance_type_dispatcher);
}