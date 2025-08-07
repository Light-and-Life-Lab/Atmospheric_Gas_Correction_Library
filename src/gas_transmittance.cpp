<%
cfg['compiler_args'] = ['-std=c++17', '-O3', '-fopenmp']
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

Transmittance_Record_PY ozone_transmittance(const L1_Record_PY& l1_rec, const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.ozone_absorption_cross_section = static_cast<double*>(l1_rec.ozone_absorption_cross_section.request().ptr);
    l1_rec_c.ozone_concentration = static_cast<double*>(l1_rec.ozone_concentration.request().ptr);
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

    ozone_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY co2_transmittance(const L1_Record_PY& l1_rec, const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.co2_transmittance = static_cast<double*>(l1_rec.co2_transmittance.request().ptr);
    l1_rec_c.air_mass_factor_mixed_gases = static_cast<double*>(l1_rec.air_mass_factor_mixed_gases.request().ptr);
    l1_rec_c.num_amf_grid_points = l1_rec.num_amf_grid_points;
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

    co2_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY co_transmittance(const L1_Record_PY& l1_rec, const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.co_transmittance = static_cast<double*>(l1_rec.co_transmittance.request().ptr);
    l1_rec_c.air_mass_factor_mixed_gases = static_cast<double*>(l1_rec.air_mass_factor_mixed_gases.request().ptr);
    l1_rec_c.num_amf_grid_points = l1_rec.num_amf_grid_points;
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

    co_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY ch4_transmittance(const L1_Record_PY& l1_rec, const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.ch4_transmittance = static_cast<double*>(l1_rec.ch4_transmittance.request().ptr);
    l1_rec_c.air_mass_factor_mixed_gases = static_cast<double*>(l1_rec.air_mass_factor_mixed_gases.request().ptr);
    l1_rec_c.num_amf_grid_points = l1_rec.num_amf_grid_points;
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

    ch4_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY no2_transmittance(const L1_Record_PY& l1_rec, const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};

    l1_rec_c.no2_absorption_cross_section = static_cast<double*>(l1_rec.no2_absorption_cross_section.request().ptr);
    l1_rec_c.fraction_tropospheric_no2_above_200m = static_cast<double*>(l1_rec.fraction_tropospheric_no2_above_200m.request().ptr);
    l1_rec_c.tropospheric_no2_concentration = static_cast<double*>(l1_rec.tropospheric_no2_concentration.request().ptr);
    l1_rec_c.stratospheric_no2_concentration = static_cast<double*>(l1_rec.stratospheric_no2_concentration.request().ptr);
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

    no2_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


PYBIND11_MODULE(gas_transmittance, m) 
{
    py::class_<L1_Record_PY>(m, "L1_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("ozone_absorption_cross_section", &L1_Record_PY::ozone_absorption_cross_section)
        .def_readwrite("ozone_concentration", &L1_Record_PY::ozone_concentration)
        .def_readwrite("no2_absorption_cross_section", &L1_Record_PY::no2_absorption_cross_section)
        .def_readwrite("fraction_tropospheric_no2_above_200m", &L1_Record_PY::fraction_tropospheric_no2_above_200m)
        .def_readwrite("tropospheric_no2_concentration", &L1_Record_PY::tropospheric_no2_concentration)
        .def_readwrite("stratospheric_no2_concentration", &L1_Record_PY::stratospheric_no2_concentration)
        .def_readwrite("co2_transmittance", &L1_Record_PY::co2_transmittance)
        .def_readwrite("co_transmittance", &L1_Record_PY::co_transmittance)
        .def_readwrite("ch4_transmittance", &L1_Record_PY::ch4_transmittance)
        .def_readwrite("air_mass_factor_mixed_gases", &L1_Record_PY::air_mass_factor_mixed_gases)
        .def_readwrite("num_amf_grid_points", &L1_Record_PY::num_amf_grid_points)
        .def_readwrite("cos_solar_zenith", &L1_Record_PY::cos_solar_zenith)
        .def_readwrite("cos_sensor_zenith", &L1_Record_PY::cos_sensor_zenith)
        .def_readwrite("num_pixels", &L1_Record_PY::num_pixels)
        .def_readwrite("num_wavelengths", &L1_Record_PY::num_wavelengths);

    py::class_<Transmittance_Record_PY>(m, "Transmittance_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("gas_transmittance_solar_zenith", &Transmittance_Record_PY::gas_transmittance_solar_zenith)
        .def_readwrite("gas_transmittance_sensor_zenith", &Transmittance_Record_PY::gas_transmittance_sensor_zenith)
        .def_readwrite("gas_transmittance_total", &Transmittance_Record_PY::gas_transmittance_total);

    m.def("ozone_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&ozone_transmittance));
    m.def("co2_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&co2_transmittance));
    m.def("co_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&co_transmittance));
    m.def("ch4_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&ch4_transmittance));
    m.def("no2_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&no2_transmittance));
}