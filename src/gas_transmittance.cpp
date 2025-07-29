<%
cfg['compiler_args'] = ['-std=c++17', '-O3', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
setup_pybind11(cfg)
%>

#include <string>
#include <iostream>
#include <cmath>
#include <tuple>
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "gas_transmittance.hpp"
#include "L1_record.h"

namespace py = pybind11;

std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> ozone_transmittance(
    const py::array_t<double, py::array::c_style>& k_oz,
    const py::array_t<double, py::array::c_style>& l1b_oz,
    const py::array_t<double, py::array::c_style>& l1b_csolz,
    const py::array_t<double, py::array::c_style>& l1b_csenz,
    const bool do_amf_correction) 
{
    double* k_oz_ptr = static_cast<double*>(k_oz.request().ptr);
    double* l1b_oz_ptr = static_cast<double*>(l1b_oz.request().ptr);
    double* l1b_csolz_ptr = static_cast<double*>(l1b_csolz.request().ptr);
    double* l1b_csenz_ptr = static_cast<double*>(l1b_csenz.request().ptr);

    int n_rows = l1b_oz.request().shape[0];
    int num_pixels = n_rows;
    int n_cols = k_oz.request().shape[0];
    int num_wavelengths = n_cols;

    py::array_t<double> tg_sol = allocate_output_array<double>(n_rows, n_cols);
    py::array_t<double> tg_sen = allocate_output_array<double>(n_rows, n_cols);
    py::array_t<double> tg = allocate_output_array<double>(n_rows, n_cols);

    double* tg_sol_ptr = static_cast<double*>(tg_sol.request().ptr);
    double* tg_sen_ptr = static_cast<double*>(tg_sen.request().ptr);
    double* tg_ptr = static_cast<double*>(tg.request().ptr);

    // TODO: Possibly factor out the core for loop that can be written entirely in C into a separate file. 
    // May be able to offer a purely C-based API that just takes the pointers to the arrays in addition to the Python API

    // Input arrays have values for each pixel at each wavelength
    // If the array is considered as a 2D matrix then each row contains the values at every wavelength for a given pixel
    #pragma omp parallel for
    for (int ip = 0; ip < num_pixels; ip++)
    {
        int row_offset = ip*num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < num_wavelengths; iw++) 
        {
            double tau_oz = l1b_oz_ptr[ip] * k_oz_ptr[iw];
            tg_sol_ptr[row_offset + iw] = exp(-(tau_oz / l1b_csolz_ptr[ip]));

            if (do_amf_correction) 
            {
                tg_ptr[row_offset + iw] = exp(-tau_oz * (1.0/l1b_csolz_ptr[ip] + 1.0/l1b_csenz_ptr[ip]));
            } 
            else 
            {
                tg_sen_ptr[row_offset + iw] = exp(-(tau_oz / l1b_csenz_ptr[ip]));
            }
        }
    }

    return std::make_tuple(tg_sol, tg_sen, tg);
}


std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> ozone_transmittance(
    const L1_Record& l1_rec,
    const bool do_amf_correction) 
{
    return ozone_transmittance(l1_rec.k_oz, l1_rec.l1b_oz, l1_rec.l1b_csolz, l1_rec.l1b_csenz, do_amf_correction);
}


std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> no2_transmittance(
    const py::array_t<double, py::array::c_style>& k_no2, 
    const py::array_t<double, py::array::c_style>& l1b_no2_frac, 
    const py::array_t<double, py::array::c_style>& l1b_no2_tropo, 
    const py::array_t<double, py::array::c_style>& l1b_no2_strat,
    const py::array_t<double, py::array::c_style>& l1b_csolz, 
    const py::array_t<double, py::array::c_style>& l1b_csenz,
    const bool do_amf_correction) 
{
    double* k_no2_ptr = static_cast<double*>(k_no2.request().ptr);
    double* l1b_no2_frac_ptr = static_cast<double*>(l1b_no2_frac.request().ptr);
    double* l1b_no2_tropo_ptr = static_cast<double*>(l1b_no2_tropo.request().ptr);
    double* l1b_no2_strat_ptr = static_cast<double*>(l1b_no2_strat.request().ptr);
    double* l1b_csolz_ptr = static_cast<double*>(l1b_csolz.request().ptr);
    double* l1b_csenz_ptr = static_cast<double*>(l1b_csenz.request().ptr);

    int n_rows = l1b_no2_tropo.request().shape[0];
    int num_pixels = n_rows;
    int n_cols = k_no2.request().shape[0];
    int num_wavelengths = n_cols;

    py::array_t<double> tg_sol = allocate_output_array<double>(n_rows, n_cols);
    py::array_t<double> tg_sen = allocate_output_array<double>(n_rows, n_cols);
    py::array_t<double> tg = allocate_output_array<double>(n_rows, n_cols);

    double* tg_sol_ptr = static_cast<double*>(tg_sol.request().ptr);
    double* tg_sen_ptr = static_cast<double*>(tg_sen.request().ptr);
    double* tg_ptr = static_cast<double*>(tg.request().ptr);

    #pragma omp parallel for
    for (int ip = 0; ip < num_pixels; ip++)
    {
        double sec0 = 1.0 / l1b_csolz_ptr[ip];
        double sec = 1.0 / l1b_csenz_ptr[ip];
        double no2_tr200{0.0};

        if (l1b_no2_tropo_ptr[ip] > 0.0)
        {
            /* compute tropo no2 above 200m (Z.Ahmad)
            no2_tr200 = exp(12.6615 + 0.61676*log(no2_tropo));
            new, location-dependent method */
            no2_tr200 = l1b_no2_frac_ptr[ip] * l1b_no2_tropo_ptr[ip];
        }

        int row_offset = ip*num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < num_wavelengths; iw++) 
        {
            if (k_no2_ptr[iw] > 0.0) 
            {
                double a_285 = k_no2_ptr[iw] * (1.0 - 0.003 * (285.0 - 294.0));
                double a_225 = k_no2_ptr[iw] * (1.0 - 0.003 * (225.0 - 294.0));

                double tau_to200 = a_285 * no2_tr200 + a_225 * l1b_no2_strat_ptr[ip];

                tg_sol_ptr[row_offset + iw] = exp(-(tau_to200 * sec0));

                if (do_amf_correction) 
                {
                    tg_ptr[row_offset + iw] = exp(-(tau_to200 * (sec + sec0)));
                }
                else
                {
                    tg_sen_ptr[row_offset + iw] = exp(-(tau_to200 * sec));
                }
            }
        }
    }

    return std::make_tuple(tg_sol, tg_sen, tg);
}


std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> no2_transmittance(
    const L1_Record& l1_rec,
    const bool do_amf_correction) 
{
    return no2_transmittance(
        l1_rec.k_no2, 
        l1_rec.l1b_no2_frac,
        l1_rec.l1b_no2_tropo,
        l1_rec.l1b_no2_strat,
        l1_rec.l1b_csolz,
        l1_rec.l1b_csenz,
        do_amf_correction
        );
}


PYBIND11_MODULE(gas_transmittance, m) 
{
    py::class_<L1_Record>(m, "L1_Record", py::module_local())
        .def(py::init<>())
        .def_readwrite("k_oz", &L1_Record::k_oz)
        .def_readwrite("l1b_oz", &L1_Record::l1b_oz)
        .def_readwrite("k_no2", &L1_Record::k_no2)
        .def_readwrite("l1b_no2_frac", &L1_Record::l1b_no2_frac)
        .def_readwrite("l1b_no2_tropo", &L1_Record::l1b_no2_tropo)
        .def_readwrite("l1b_no2_strat", &L1_Record::l1b_no2_strat)
        .def_readwrite("l1b_csolz", &L1_Record::l1b_csolz)
        .def_readwrite("l1b_csenz", &L1_Record::l1b_csenz);

    m.def("ozone_transmittance", py::overload_cast<const L1_Record&, bool>(&ozone_transmittance));
    m.def("ozone_transmittance", py::overload_cast<
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&, 
        bool> (&ozone_transmittance));

    m.def("no2_transmittance", py::overload_cast<const L1_Record&, bool>(&no2_transmittance));
    m.def("no2_transmittance", py::overload_cast<
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        const py::array_t<double, py::array::c_style>&,
        bool>(&no2_transmittance));
}