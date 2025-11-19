import numpy as np
import matplotlib.pyplot as plt
import h5py
from pyhdf.V import *
from pyhdf.HDF import *
from pyhdf.SD import *
from scipy import interpolate
import time
import pytest

import bin.gas_transmittance as gas_transmittance

save_transmittances = False

@pytest.fixture
def read_PACE_geometry_data():
    csolz = np.load('test/PACE/csolz.npy')
    csenz = np.load('test/PACE/csenz.npy')

    return csolz, csenz

@pytest.fixture
def read_PACE_data():
    l1_data = gas_transmittance.L1_Data()
    start_line = 0
    start_pixel = 0
    end_line = 100
    end_pixel = 100

    with h5py.File('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', 'r') as f:
        solar_zenith = 0.01*np.array(f['/geolocation_data/solar_zenith'][start_line:end_line, start_pixel:end_pixel])
        sensor_zenith = 0.01*np.array(f['/geolocation_data/sensor_zenith'][start_line:end_line, start_pixel:end_pixel])
        l1_data.cos_solar_zenith = np.cos(np.deg2rad(solar_zenith))
        l1_data.cos_sensor_zenith = np.cos(np.deg2rad(sensor_zenith))
        blue_wavelengths = np.array(f['/sensor_band_parameters/blue_wavelength'][1:])
        red_wavelengths = np.array(f['/sensor_band_parameters/red_wavelength'][3:])
        sensor_wavelengths = np.zeros(len(blue_wavelengths) + len(red_wavelengths))
        sensor_wavelengths[0:len(blue_wavelengths)] = blue_wavelengths
        sensor_wavelengths[len(blue_wavelengths):] = red_wavelengths

        blue_rhot = np.array(f['/observation_data/rhot_blue'][1:, start_line:end_line, start_pixel:end_pixel])
        red_rhot = np.array(f['/observation_data/rhot_red'][3:, start_line:end_line, start_pixel:end_pixel])
        assert(blue_rhot.shape[1] == red_rhot.shape[1])
        assert(blue_rhot.shape[2] == red_rhot.shape[2])

        rhot = np.zeros((blue_rhot.shape[0] + red_rhot.shape[0], blue_rhot.shape[1], blue_rhot.shape[2]))
        rhot[0:len(blue_wavelengths), :, :] = blue_rhot
        rhot[len(blue_wavelengths):, :, :] = red_rhot
        rhot = np.rollaxis(rhot, 0, 3)
        l1_data.reflectance = rhot/np.pi*l1_data.cos_sensor_zenith[:, :, None]

        assert(len(sensor_wavelengths) == l1_data.reflectance.shape[2])

        l1_data.wavelengths = sensor_wavelengths
        l1_data.num_pixels = len(l1_data.cos_solar_zenith.flatten())
        l1_data.num_wavelengths = len(l1_data.wavelengths)

    return l1_data


@pytest.fixture
def read_OCSSW_lat_lon():
    with h5py.File('test/PACE/ozone/PACE_OCI.20240411T182012.L2.ozone.nc', 'r') as f:
        OCSSW_lat = np.array(f['/navigation_data/latitude'])
        OCSSW_lon = np.array(f['/navigation_data/longitude'])

    return OCSSW_lat, OCSSW_lon


@pytest.fixture
def read_OCSMART_lat_lon():
    OCSMART_lat = np.load('test/PACE/l1b_lat.npy')
    OCSMART_lon = np.load('test/PACE/l1b_lon.npy')

    return OCSMART_lat, OCSMART_lon


@pytest.fixture
def read_gas_transmittance_table():
    gas_transmittance_table = gas_transmittance.Gas_Transmittance_Lookup_Table()

    netcdf_variables = {
            'air_mass_factor_mixed' : "air_mass_factor_mixed_gases",
            'air_mass_factor_wv' : "air_mass_factor_water_vapor",
            'carbon_dioxide_transmittance' : "co2_transmittance",
            'carbon_monoxide_transmittance' : "co_transmittance",
            'methane_transmittance' : "ch4_transmittance",
            'nitrous_oxide_transmittance' : "n2o_transmittance",
            'oxygen_transmittance' : "o2_transmittance",
            'water_vapor' : "water_vapor_concentration",
            'water_vapor_transmittance' : "h2o_transmittance",
            'wavelength' : "wavelengths",
        }
    
    netcdf_dimensions = {
        'n_air_mass_factor' : "num_amf_grid_points",
        'n_water_vapor' : "num_water_vapor_concentrations",
        'nmodels' : 'num_models',
        'nwavelengths' : 'num_wavelengths',
    }

    model_map = {
        "Tropical": 0,
        "MidLatSummer" : 1,
        "MidLatWinter" : 2,
        "SubarcticSummer" : 3,
        "SubarcticWinter" : 4,
        "USstandard62" : 5
    }
        
    with h5py.File('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc', 'r') as f:
        for netcdf_node, var_name in netcdf_variables.items():
            if netcdf_node in f:
                var_value = np.array(f[netcdf_node])
                setattr(gas_transmittance_table, var_name, var_value)

        for netcdf_node, var_name in netcdf_dimensions.items():
            if netcdf_node in f:
                var_value = np.array(f[netcdf_node])
                setattr(gas_transmittance_table, var_name, len(var_value))

        setattr(gas_transmittance_table, 'model', model_map["USstandard62"])

    return gas_transmittance_table


@pytest.fixture
def read_ozone_ancillary_data():
    koz = np.load('test/PACE/ozone/koz.npy')
    ozone_concentration = np.load('test/PACE/ozone/oz_concentration.npy')

    return koz, ozone_concentration


@pytest.fixture
def read_OCSMART_ozone_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/PACE/ozone/tg_sol_oz.npy')
    tg_sen_ocsmart = np.load('test/PACE/ozone/tg_sen_oz.npy')
    sensor_wavelengths = np.load('test/PACE/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture
def read_OCSSW_ozone_transmittance_benchmark_data():
    with h5py.File('test/PACE/ozone/PACE_OCI.20240411T182012.L2.ozone.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_no2_ancillary_data():
    kno2 = np.load('test/PACE/no2/k_no2.npy')
    no2_frac = np.load('test/PACE/no2/no2_frac.npy')
    no2_tropo = np.load('test/PACE/no2/no2_tropo.npy')
    no2_strat = np.load('test/PACE/no2/no2_strat.npy')

    return kno2, no2_frac, no2_tropo, no2_strat


@pytest.fixture
def read_OCSMART_no2_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/PACE/no2/tg_sol_no2.npy')
    tg_sen_ocsmart = np.load('test/PACE/no2/tg_sen_no2.npy')
    sensor_wavelengths = np.load('test/PACE/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture
def read_OCSSW_no2_transmittance_benchmark_data():
    with h5py.File('test/PACE/no2/PACE_OCI.20240411T182012.L2.no2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_co2_transmittance_benchmark_data():
    with h5py.File('test/PACE/co2/PACE_OCI.20240411T182012.L2.co2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_co_transmittance_benchmark_data():
    with h5py.File('test/PACE/co/PACE_OCI.20240411T182012.L2.co.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_ch4_transmittance_benchmark_data():
    with h5py.File('test/PACE/ch4/PACE_OCI.20240411T182012.L2.ch4.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_n2o_transmittance_benchmark_data():
    with h5py.File('test/PACE/n2o/PACE_OCI.20240411T182012.L2.n2o.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_o2_transmittance_benchmark_data():
    with h5py.File('test/PACE/o2/PACE_OCI.20240411T182012.L2.o2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture
def read_OCSSW_h2o_transmittance_benchmark_data():
    with h5py.File('test/PACE/h2o/PACE_OCI.20240411T182012.L2.h2o.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


def save_gas_transmittances(gas_transmittances, gas_subfolder):
    output_path = 'test/PACE/' + gas_subfolder + '/'
    np.save(output_path + 'solar_zenith.npy', gas_transmittances.solar_zenith)
    np.save(output_path + 'sensor_zenith.npy', gas_transmittances.sensor_zenith)
    np.save(output_path + 'total.npy', gas_transmittances.total)


def load_gas_transmittances(gas_subfolder):
    output_path = 'test/PACE/' + gas_subfolder + '/'
    solar_zenith_saved = np.load(output_path + 'solar_zenith.npy')
    sensor_zenith_saved = np.load(output_path + 'sensor_zenith.npy')
    total_saved = np.load(output_path + 'total.npy')

    return solar_zenith_saved, sensor_zenith_saved, total_saved


# @pytest.mark.skip()
def test_ozone_OCSSW(read_ozone_ancillary_data, 
                     read_PACE_geometry_data, 
                     read_OCSMART_ozone_transmittance_benchmark_data, 
                     read_OCSSW_ozone_transmittance_benchmark_data):
    ancillary_data = gas_transmittance.Ancillary_Data()
    ancillary_data.ozone_absorption_cross_section, ancillary_data.ozone_concentration = read_ozone_ancillary_data

    # TODO: Adapt approach from OCSMART ancillary.py to read from ancillary files directly instead of .npy files
    # TODO: This will require interpolation of ozone data to l1b grid
    l1_data = gas_transmittance.L1_Data()
    l1_data.cos_solar_zenith, l1_data.cos_sensor_zenith = read_PACE_geometry_data
    l1_data.num_pixels = l1_data.cos_solar_zenith.shape[0] * l1_data.cos_solar_zenith.shape[1]
    l1_data.num_wavelengths = len(ancillary_data.ozone_absorption_cross_section)
    do_amf_correction = False

    gas_transmittances = gas_transmittance.ozone_transmittance(l1_data, ancillary_data, do_amf_correction)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_ozone_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((1710, 1272, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((1710, 1272, 197))

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_ozone_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'ozone')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('ozone')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths, tg_sen_ocsmart[0, 0, :], '-g')
    plt.plot(sensor_wavelengths, tg_sol_ocsmart[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths, tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith'])
    plt.savefig('test/PACE/ozone/transmittance_comparison.png')


# @pytest.mark.skip()
def test_co2_OCSSW(read_gas_transmittance_table, 
                    read_PACE_data,  
                    read_OCSSW_co2_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = True

    gas_transmittances = gas_transmittance.co2_transmittance(l1_data, gas_transmittance_table, do_amf_correction)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co2_transmittance_benchmark_data

    sensor_wavelengths = l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'co2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('co2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/co2/transmittance_comparison.png')


# @pytest.mark.skip()
def test_co_OCSSW(read_gas_transmittance_table, 
                    read_PACE_data,  
                    read_OCSSW_co_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = False

    gas_transmittances = gas_transmittance.co_transmittance(l1_data, gas_transmittance_table, do_amf_correction)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co_transmittance_benchmark_data

    sensor_wavelengths = l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'co')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('co')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/co/transmittance_comparison.png')


# @pytest.mark.skip()
def test_ch4_OCSSW(read_gas_transmittance_table, 
                    read_PACE_data,  
                    read_OCSSW_ch4_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = False

    gas_transmittances = gas_transmittance.ch4_transmittance(l1_data, gas_transmittance_table, do_amf_correction)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_ch4_transmittance_benchmark_data

    sensor_wavelengths = l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'ch4')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('ch4')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/ch4/transmittance_comparison.png')


# @pytest.mark.skip()
def test_n2o_OCSSW(read_gas_transmittance_table, 
                    read_PACE_data,  
                    read_OCSSW_n2o_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = False

    gas_transmittances = gas_transmittance.n2o_transmittance(l1_data, gas_transmittance_table, do_amf_correction)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_n2o_transmittance_benchmark_data

    sensor_wavelengths = l1_data.wavelengths

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'n2o')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('n2o')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/n2o/transmittance_comparison.png')


# @pytest.mark.skip()
def test_no2_OCSSW(read_no2_ancillary_data, 
                     read_PACE_geometry_data, 
                     read_OCSMART_no2_transmittance_benchmark_data, 
                     read_OCSSW_no2_transmittance_benchmark_data):
    ancillary_data = gas_transmittance.Ancillary_Data()

    ancillary_data.no2_absorption_cross_section, \
    ancillary_data.fraction_tropospheric_no2_above_200m, \
    ancillary_data.tropospheric_no2_concentration, \
    ancillary_data.stratospheric_no2_concentration = read_no2_ancillary_data

    # TODO: Adapt approach from OCSMART ancillary.py to read from ancillary files directly instead of .npy files
    # TODO: This will require interpolation of ozone data to l1b grid
    l1_data = gas_transmittance.L1_Data()
    l1_data.cos_solar_zenith, l1_data.cos_sensor_zenith = read_PACE_geometry_data
    l1_data.num_pixels = l1_data.cos_solar_zenith.shape[0] * l1_data.cos_solar_zenith.shape[1]
    l1_data.num_wavelengths = len(ancillary_data.no2_absorption_cross_section)
    do_amf_correction = False

    gas_transmittances = gas_transmittance.no2_transmittance(l1_data, ancillary_data, do_amf_correction)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_no2_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((1710, 1272, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((1710, 1272, 197))

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_no2_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'no2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('no2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(sensor_wavelengths, tg_sen_ocsmart[0, 0, :], '-g')
    plt.plot(sensor_wavelengths, tg_sol_ocsmart[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(sensor_wavelengths, tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith'])
    plt.savefig('test/PACE/no2/transmittance_comparison.png')


# @pytest.mark.skip()
def test_o2_OCSSW(read_gas_transmittance_table,
                  read_PACE_data,
                  read_OCSSW_o2_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = True

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.o2_transmittance, axis = 0)
    o2_transmittance_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.o2_transmittance = o2_transmittance_sensor_wavelengths

    oxygen_A_band_option = gas_transmittance.Oxygen_A_Band_Option.TRANSMITTANCE_TABLE

    start_time = time.perf_counter()
    gas_transmittances = gas_transmittance.o2_transmittance(l1_data, gas_transmittance_table, do_amf_correction, oxygen_A_band_option)
    end_time = time.perf_counter()

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_o2_transmittance_benchmark_data
    
    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'o2')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('o2')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([750, 800])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/o2/transmittance_comparison_oxaband_opt_2.png')


# @pytest.mark.skip()
def test_h2o_OCSSW(read_gas_transmittance_table,
                   read_PACE_data,
                   read_no2_ancillary_data,
                   read_OCSSW_h2o_transmittance_benchmark_data):
    gas_transmittance_table = read_gas_transmittance_table

    l1_data = read_PACE_data
    do_amf_correction = True

    ancillary_data = gas_transmittance.Ancillary_Data()

    ancillary_data.no2_absorption_cross_section, \
    ancillary_data.fraction_tropospheric_no2_above_200m, \
    ancillary_data.tropospheric_no2_concentration, \
    ancillary_data.stratospheric_no2_concentration = read_no2_ancillary_data

    ancillary_data.precipitable_water = np.zeros(l1_data.cos_solar_zenith.size)
    ancillary_data.water_vapor_bands = np.array([782, 817, 857], dtype=np.float64)
    ancillary_data.num_water_vapor_bands = ancillary_data.water_vapor_bands.size

    f = interpolate.interp1d(gas_transmittance_table.wavelengths, gas_transmittance_table.h2o_transmittance, axis = 1)
    h2o_transmittance_at_sensor_wavelengths = f(l1_data.wavelengths)

    gas_transmittance_table.num_wavelengths = len(l1_data.wavelengths)
    gas_transmittance_table.h2o_transmittance = h2o_transmittance_at_sensor_wavelengths

    use_gas_transmittance_table = True

    start_time = time.perf_counter()
    gas_transmittances = gas_transmittance.h2o_transmittance(l1_data, ancillary_data, gas_transmittance_table, do_amf_correction, use_gas_transmittance_table)
    end_time = time.perf_counter()

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_h2o_transmittance_benchmark_data

    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'h2o')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('h2o')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d[:], tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d[:], tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([500, 900])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/h2o/transmittance_comparison.png')