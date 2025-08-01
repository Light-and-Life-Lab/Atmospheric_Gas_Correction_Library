#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

struct L1_Record_PY
{
public:
    L1_Record_PY(){};
    
    pybind11::array_t<double, pybind11::array::c_style> k_oz{};
    pybind11::array_t<double, pybind11::array::c_style> l1b_oz{}; 
    pybind11::array_t<double, pybind11::array::c_style> k_no2{};
    pybind11::array_t<double, pybind11::array::c_style> l1b_no2_frac{};
    pybind11::array_t<double, pybind11::array::c_style> l1b_no2_tropo{};
    pybind11::array_t<double, pybind11::array::c_style> l1b_no2_strat{};
    pybind11::array_t<double, pybind11::array::c_style> l1b_csolz{}; 
    pybind11::array_t<double, pybind11::array::c_style> l1b_csenz{};
};