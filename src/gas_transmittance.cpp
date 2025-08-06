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

Transmittance_Record_PY ozone_transmittance(
    const L1_Record_PY& l1_rec,
    const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.k_oz = static_cast<double*>(l1_rec.k_oz.request().ptr);
    l1_rec_c.l1b_oz = static_cast<double*>(l1_rec.l1b_oz.request().ptr);
    l1_rec_c.l1b_csolz = static_cast<double*>(l1_rec.l1b_csolz.request().ptr);
    l1_rec_c.l1b_csenz = static_cast<double*>(l1_rec.l1b_csenz.request().ptr);

    int n_rows = l1_rec.l1b_oz.request().shape[0];
    l1_rec_c.num_pixels = n_rows;
    int n_cols = l1_rec.k_oz.request().shape[0];
    l1_rec_c.num_wavelengths = n_cols;

    Transmittance_Record_PY t_rec{};

    t_rec.tg_sol = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg_sen = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.tg_sol = static_cast<double*>(t_rec.tg_sol.request().ptr);
    t_rec_c.tg_sen = static_cast<double*>(t_rec.tg_sen.request().ptr);
    t_rec_c.tg = static_cast<double*>(t_rec.tg.request().ptr);

    ozone_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY co_transmittance(
    const L1_Record_PY& l1_rec,
    const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};
    
    l1_rec_c.t_co = static_cast<double*>(l1_rec.t_co.request().ptr);
    l1_rec_c.amf_mixed = static_cast<double*>(l1_rec.amf_mixed.request().ptr);
    l1_rec_c.num_airmass = l1_rec.num_airmass;
    l1_rec_c.l1b_csolz = static_cast<double*>(l1_rec.l1b_csolz.request().ptr);
    l1_rec_c.l1b_csenz = static_cast<double*>(l1_rec.l1b_csenz.request().ptr);

    int n_rows = l1_rec.l1b_csolz.request().shape[0];
    l1_rec_c.num_pixels = n_rows;
    int n_cols = l1_rec.t_co.request().shape[0];
    l1_rec_c.num_wavelengths = n_cols;

    Transmittance_Record_PY t_rec{};

    t_rec.tg_sol = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg_sen = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.tg_sol = static_cast<double*>(t_rec.tg_sol.request().ptr);
    t_rec_c.tg_sen = static_cast<double*>(t_rec.tg_sen.request().ptr);
    t_rec_c.tg = static_cast<double*>(t_rec.tg.request().ptr);

    co_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


Transmittance_Record_PY no2_transmittance(
    const L1_Record_PY& l1_rec,
    const bool do_amf_correction) 
{
    L1_Record l1_rec_c{};

    l1_rec_c.k_no2 = static_cast<double*>(l1_rec.k_no2.request().ptr);
    l1_rec_c.l1b_no2_frac = static_cast<double*>(l1_rec.l1b_no2_frac.request().ptr);
    l1_rec_c.l1b_no2_tropo = static_cast<double*>(l1_rec.l1b_no2_tropo.request().ptr);
    l1_rec_c.l1b_no2_strat = static_cast<double*>(l1_rec.l1b_no2_strat.request().ptr);
    l1_rec_c.l1b_csolz = static_cast<double*>(l1_rec.l1b_csolz.request().ptr);
    l1_rec_c.l1b_csenz = static_cast<double*>(l1_rec.l1b_csenz.request().ptr);

    int n_rows = l1_rec.l1b_no2_tropo.request().shape[0];
    l1_rec_c.num_pixels = n_rows;
    int n_cols = l1_rec.k_no2.request().shape[0];
    l1_rec_c.num_wavelengths = n_cols;

    Transmittance_Record_PY t_rec{};

    t_rec.tg_sol = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg_sen = allocate_output_array<double>(n_rows, n_cols);
    t_rec.tg = allocate_output_array<double>(n_rows, n_cols);

    Transmittance_Record t_rec_c{};

    t_rec_c.tg_sol = static_cast<double*>(t_rec.tg_sol.request().ptr);
    t_rec_c.tg_sen = static_cast<double*>(t_rec.tg_sen.request().ptr);
    t_rec_c.tg = static_cast<double*>(t_rec.tg.request().ptr);

    no2_transmittance(&l1_rec_c, &t_rec_c, do_amf_correction);

    return t_rec;
}


PYBIND11_MODULE(gas_transmittance, m) 
{
    py::class_<L1_Record_PY>(m, "L1_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("k_oz", &L1_Record_PY::k_oz)
        .def_readwrite("l1b_oz", &L1_Record_PY::l1b_oz)
        .def_readwrite("k_no2", &L1_Record_PY::k_no2)
        .def_readwrite("l1b_no2_frac", &L1_Record_PY::l1b_no2_frac)
        .def_readwrite("l1b_no2_tropo", &L1_Record_PY::l1b_no2_tropo)
        .def_readwrite("l1b_no2_strat", &L1_Record_PY::l1b_no2_strat)
        .def_readwrite("t_co", &L1_Record_PY::t_co)
        .def_readwrite("amf_mixed", &L1_Record_PY::amf_mixed)
        .def_readwrite("num_airmass", &L1_Record_PY::num_airmass)
        .def_readwrite("l1b_csolz", &L1_Record_PY::l1b_csolz)
        .def_readwrite("l1b_csenz", &L1_Record_PY::l1b_csenz);

    py::class_<Transmittance_Record_PY>(m, "Transmittance_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("tg_sol", &Transmittance_Record_PY::tg_sol)
        .def_readwrite("tg_sen", &Transmittance_Record_PY::tg_sen)
        .def_readwrite("tg", &Transmittance_Record_PY::tg);

    m.def("ozone_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&ozone_transmittance));
    m.def("co_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&co_transmittance));
    m.def("no2_transmittance", py::overload_cast<const L1_Record_PY&, bool>(&no2_transmittance));
}