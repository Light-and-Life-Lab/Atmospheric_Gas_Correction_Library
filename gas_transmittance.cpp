<%
cfg['compiler_args'] = ['-std=c++17']
setup_pybind11(cfg)
%>

#include <string>
#include <iostream>
#include <cmath>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "gas_transmittance.hpp"

namespace py = pybind11;

std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>> ozone_transmittance(py::array_t<double> k_oz, 
                         py::array_t<double> l1b_oz, 
                         py::array_t<double> l1b_solz, 
                         py::array_t<double> l1b_senz,
                         bool do_amf_correction) 
{
    double* k_oz_ptr = static_cast<double*>(k_oz.request().ptr);
    double* l1b_oz_ptr = static_cast<double*>(l1b_oz.request().ptr);
    double* l1b_solz_ptr = static_cast<double*>(l1b_solz.request().ptr);
    double* l1b_senz_ptr = static_cast<double*>(l1b_senz.request().ptr);

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

    double tau_oz{std::numeric_limits<double>::quiet_NaN()};

    // Input arrays have values for each pixel at each wavelength
    // If the array is considered as a 2D matrix then each row contains the values at every wavelength for a given pixel
    py::print("Computing transmittances");
    for (int ip = 0; ip < num_pixels; ip++)
    {
        py::print("\nip = ", ip);
        int row_offset = ip*num_wavelengths; // Each row represents a single pixel and has num_wavelengths elements

        for (int iw = 0; iw < num_wavelengths; iw++) 
        {
            py::print("iw = ", iw);
            tau_oz = l1b_oz_ptr[ip] * k_oz_ptr[iw];
            tg_sol_ptr[row_offset + iw] = exp(-(tau_oz / l1b_solz_ptr[ip]));

            if (do_amf_correction) 
            {
                tg_ptr[row_offset + iw] = exp(-tau_oz * (1.0/l1b_solz_ptr[ip] + 1.0/l1b_senz_ptr[ip]));
            } 
            else 
            {
                tg_sen_ptr[row_offset + iw] = exp(-(tau_oz / l1b_senz_ptr[ip]));
            }
        }
    }

    return std::make_tuple(tg_sol, tg_sen, tg);
}

PYBIND11_MODULE(gas_transmittance, m) 
{
    m.def("ozone_transmittance", &ozone_transmittance);
}