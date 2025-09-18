<%
cfg['compiler_args'] = ['-std=c++2a', '-O3', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
cfg['sources'] = ['gas_transmittance_impl.cpp']
setup_pybind11(cfg)
%>

#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "allocate_output_array.hpp"
#include "pybind_interface_types.hpp"
#include "gas_transmittance.h"

namespace py = pybind11;

Transmittance_Record_PY ozone_transmittance(const L1_Record_PY& l1_rec, const Ancillary_Data_PY& ancillary_data, const bool do_amf_correction) 
{
    Ancillary_Data ancillary_data_c{};
    
    ancillary_data_c.ozone_absorption_cross_section = static_cast<double*>(ancillary_data.ozone_absorption_cross_section.request().ptr);
    ancillary_data_c.ozone_concentration = static_cast<double*>(ancillary_data.ozone_concentration.request().ptr);
    
    L1_Record l1_rec_c{};

    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    ozone_transmittance(&l1_rec_c, &ancillary_data_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY co2_transmittance(const L1_Record_PY& l1_rec, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction) 
{
    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.co2_transmittance = static_cast<double*>(amf_table.co2_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_mixed_gases.request().ptr);
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;

    L1_Record l1_rec_c{};

    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    co2_transmittance(&l1_rec_c, &amf_table_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY co_transmittance(const L1_Record_PY& l1_rec, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction) 
{
    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.co_transmittance = static_cast<double*>(amf_table.co_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_mixed_gases.request().ptr);
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;
    
    L1_Record l1_rec_c{};

    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    co_transmittance(&l1_rec_c, &amf_table_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY ch4_transmittance(const L1_Record_PY& l1_rec, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction) 
{
    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.ch4_transmittance = static_cast<double*>(amf_table.ch4_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_mixed_gases.request().ptr);
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;

    L1_Record l1_rec_c{};

    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    ch4_transmittance(&l1_rec_c, &amf_table_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY o2_transmittance(const L1_Record_PY& l1_rec, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction, Oxygen_A_Band_Option oxygen_A_band_option) 
{
    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.o2_transmittance = static_cast<double*>(amf_table.o2_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_mixed_gases.request().ptr);
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;
    
    L1_Record l1_rec_c{};

    l1_rec_c.Lt = static_cast<double*>(l1_rec.Lt.request().ptr);
    l1_rec_c.F0 = static_cast<double*>(l1_rec.F0.request().ptr);
    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;
    l1_rec_c.wavelengths = static_cast<double*>(l1_rec.wavelengths.request().ptr);

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    o2_transmittance(&l1_rec_c, &amf_table_c, &t_rec_c, do_amf_correction, oxygen_A_band_option);

    return t_rec;
}


Transmittance_Record_PY n2o_transmittance(const L1_Record_PY& l1_rec, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction) 
{
    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.n2o_transmittance = static_cast<double*>(amf_table.n2o_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_mixed_gases.request().ptr);
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;
    
    L1_Record l1_rec_c{};

    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    n2o_transmittance(&l1_rec_c, &amf_table_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY no2_transmittance(const L1_Record_PY& l1_rec, const Ancillary_Data_PY& ancillary_data, const bool do_amf_correction) 
{
    Ancillary_Data ancillary_data_c{};

    ancillary_data_c.no2_absorption_cross_section = static_cast<double*>(ancillary_data.no2_absorption_cross_section.request().ptr);
    ancillary_data_c.fraction_tropospheric_no2_above_200m = static_cast<double*>(ancillary_data.fraction_tropospheric_no2_above_200m.request().ptr);
    ancillary_data_c.tropospheric_no2_concentration = static_cast<double*>(ancillary_data.tropospheric_no2_concentration.request().ptr);
    ancillary_data_c.stratospheric_no2_concentration = static_cast<double*>(ancillary_data.stratospheric_no2_concentration.request().ptr);
    
    L1_Record l1_rec_c{};
    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    no2_transmittance(&l1_rec_c, &ancillary_data_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY h2o_transmittance(const L1_Record_PY& l1_rec, const Ancillary_Data_PY& ancillary_data, const Air_Mass_Factor_Lookup_Table_PY& amf_table, const bool do_amf_correction, const bool use_gas_transmittance_table) 
{
    Ancillary_Data ancillary_data_c{};

    ancillary_data_c.no2_absorption_cross_section = static_cast<double*>(ancillary_data.no2_absorption_cross_section.request().ptr);
    ancillary_data_c.fraction_tropospheric_no2_above_200m = static_cast<double*>(ancillary_data.fraction_tropospheric_no2_above_200m.request().ptr);
    ancillary_data_c.tropospheric_no2_concentration = static_cast<double*>(ancillary_data.tropospheric_no2_concentration.request().ptr);
    ancillary_data_c.stratospheric_no2_concentration = static_cast<double*>(ancillary_data.stratospheric_no2_concentration.request().ptr);

    Air_Mass_Factor_Lookup_Table amf_table_c{};

    amf_table_c.h2o_transmittance = static_cast<double*>(amf_table.h2o_transmittance.request().ptr);
    amf_table_c.air_mass_factor_mixed_gases = static_cast<double*>(amf_table.air_mass_factor_water_vapor.request().ptr);
    amf_table_c.gas_transmittance_table_wavelengths = static_cast<double*>(amf_table.gas_transmittance_table_wavelengths.request().ptr);
    amf_table_c.water_vapor_concentration = static_cast<double*>(amf_table.water_vapor_concentration.request().ptr);
    amf_table_c.num_models = amf_table.num_models;
    amf_table_c.num_gas_transmittance_wavelengths = amf_table.num_gas_transmittance_wavelengths;
    amf_table_c.num_amf_grid_points = amf_table.num_amf_grid_points;
    amf_table_c.num_water_vapor_concentrations = amf_table.num_water_vapor_concentrations;
    
    L1_Record l1_rec_c{};
    l1_rec_c.cos_solar_zenith = static_cast<double*>(l1_rec.cos_solar_zenith.request().ptr);
    l1_rec_c.cos_sensor_zenith = static_cast<double*>(l1_rec.cos_sensor_zenith.request().ptr);
    l1_rec_c.num_pixels = l1_rec.num_pixels;
    l1_rec_c.num_wavelengths = l1_rec.num_wavelengths;

    int n_rows = l1_rec_c.num_pixels;
    int n_cols = l1_rec_c.num_wavelengths;

    Transmittance_Record_PY t_rec{};

    t_rec.gas_transmittance_solar_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_sensor_zenith = allocate_output_array<double>(n_rows, n_cols);
    t_rec.gas_transmittance_total = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.gas_transmittance_solar_zenith = static_cast<double*>(t_rec.gas_transmittance_solar_zenith.request().ptr);
    t_rec_c.gas_transmittance_sensor_zenith = static_cast<double*>(t_rec.gas_transmittance_sensor_zenith.request().ptr);
    t_rec_c.gas_transmittance_total = static_cast<double*>(t_rec.gas_transmittance_total.request().ptr);

    h2o_transmittance(&l1_rec_c, &ancillary_data_c, &amf_table_c, &t_rec_c, do_amf_correction, use_gas_transmittance_table);

    return t_rec;
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
        .def_readwrite("a_h2o", &Ancillary_Data_PY::a_h2o)
        .def_readwrite("b_h2o", &Ancillary_Data_PY::b_h2o)
        .def_readwrite("c_h2o", &Ancillary_Data_PY::c_h2o)
        .def_readwrite("d_h2o", &Ancillary_Data_PY::d_h2o)
        .def_readwrite("e_h2o", &Ancillary_Data_PY::e_h2o)
        .def_readwrite("f_h2o", &Ancillary_Data_PY::f_h2o)
        .def_readwrite("g_h2o", &Ancillary_Data_PY::g_h2o)
        .def_readwrite("water_vapor_bands", &Ancillary_Data_PY::water_vapor_bands)
        .def_readwrite("num_water_vapor_bands", &Ancillary_Data_PY::num_water_vapor_bands);

    py::class_<Air_Mass_Factor_Lookup_Table_PY>(m, "Air_Mass_Factor_Lookup_Table", py::module_local())
        .def(py::init<>())
        .def_readwrite("co2_transmittance", &Air_Mass_Factor_Lookup_Table_PY::co2_transmittance)
        .def_readwrite("co_transmittance", &Air_Mass_Factor_Lookup_Table_PY::co_transmittance)
        .def_readwrite("ch4_transmittance", &Air_Mass_Factor_Lookup_Table_PY::ch4_transmittance)
        .def_readwrite("o2_transmittance", &Air_Mass_Factor_Lookup_Table_PY::o2_transmittance)
        .def_readwrite("n2o_transmittance", &Air_Mass_Factor_Lookup_Table_PY::n2o_transmittance)
        .def_readwrite("h2o_transmittance", &Air_Mass_Factor_Lookup_Table_PY::h2o_transmittance)
        .def_readwrite("model", &Air_Mass_Factor_Lookup_Table_PY::model)
        .def_readwrite("gas_transmittance_table_wavelengths", &Air_Mass_Factor_Lookup_Table_PY::gas_transmittance_table_wavelengths)
        .def_readwrite("air_mass_factor_mixed_gases", &Air_Mass_Factor_Lookup_Table_PY::air_mass_factor_mixed_gases)
        .def_readwrite("air_mass_factor_water_vapor", &Air_Mass_Factor_Lookup_Table_PY::air_mass_factor_water_vapor)
        .def_readwrite("water_vapor_concentration", &Air_Mass_Factor_Lookup_Table_PY::water_vapor_concentration)
        .def_readwrite("num_models", &Air_Mass_Factor_Lookup_Table_PY::num_models)
        .def_readwrite("num_gas_transmittance_wavelengths", &Air_Mass_Factor_Lookup_Table_PY::num_gas_transmittance_wavelengths)
        .def_readwrite("num_amf_grid_points", &Air_Mass_Factor_Lookup_Table_PY::num_amf_grid_points)
        .def_readwrite("num_water_vapor_concentrations", &Air_Mass_Factor_Lookup_Table_PY::num_water_vapor_concentrations);

    py::class_<L1_Record_PY>(m, "L1_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("Lt", &L1_Record_PY::Lt)
        .def_readwrite("F0", &L1_Record_PY::F0)
        .def_readwrite("cos_solar_zenith", &L1_Record_PY::cos_solar_zenith)
        .def_readwrite("cos_sensor_zenith", &L1_Record_PY::cos_sensor_zenith)
        .def_readwrite("num_pixels", &L1_Record_PY::num_pixels)
        .def_readwrite("num_wavelengths", &L1_Record_PY::num_wavelengths)
        .def_readwrite("wavelengths", &L1_Record_PY::wavelengths);

    py::class_<Transmittance_Record_PY>(m, "Transmittance_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("gas_transmittance_solar_zenith", &Transmittance_Record_PY::gas_transmittance_solar_zenith)
        .def_readwrite("gas_transmittance_sensor_zenith", &Transmittance_Record_PY::gas_transmittance_sensor_zenith)
        .def_readwrite("gas_transmittance_total", &Transmittance_Record_PY::gas_transmittance_total);

    py::enum_<Oxygen_A_Band_Option>(m, "Oxygen_A_Band_Option", py::module_local())
        .value("DING_GORDON", Oxygen_A_Band_Option::DING_GORDON)
        .value("NO_AMF_CORRECTION", Oxygen_A_Band_Option::NO_AMF_CORRECTION)
        .value("YES_AMF_CORRECTION", Oxygen_A_Band_Option::YES_AMF_CORRECTION);

    m.def("ozone_transmittance", py::overload_cast<const L1_Record_PY&, const Ancillary_Data_PY&, bool>(&ozone_transmittance));
    m.def("co2_transmittance", py::overload_cast<const L1_Record_PY&, const Air_Mass_Factor_Lookup_Table_PY&, bool>(&co2_transmittance));
    m.def("co_transmittance", py::overload_cast<const L1_Record_PY&, const Air_Mass_Factor_Lookup_Table_PY&, bool>(&co_transmittance));
    m.def("ch4_transmittance", py::overload_cast<const L1_Record_PY&, const Air_Mass_Factor_Lookup_Table_PY&, bool>(&ch4_transmittance));
    m.def("o2_transmittance", py::overload_cast<const L1_Record_PY&, const Air_Mass_Factor_Lookup_Table_PY& ,bool, Oxygen_A_Band_Option>(&o2_transmittance));
    m.def("n2o_transmittance", py::overload_cast<const L1_Record_PY&, const Air_Mass_Factor_Lookup_Table_PY&, bool>(&n2o_transmittance));
    m.def("no2_transmittance", py::overload_cast<const L1_Record_PY&, const Ancillary_Data_PY&, bool>(&no2_transmittance));
    m.def("h2o_transmittance", py::overload_cast<const L1_Record_PY&, const Ancillary_Data_PY&, const Air_Mass_Factor_Lookup_Table_PY&, bool, bool>(&h2o_transmittance));
}