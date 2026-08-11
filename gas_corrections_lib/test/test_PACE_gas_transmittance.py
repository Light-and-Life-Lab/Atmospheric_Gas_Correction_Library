import numpy as np
import matplotlib.pyplot as plt
import h5py
from pyhdf.V import *
from pyhdf.HDF import *
from pyhdf.SD import *
from netCDF4 import Dataset
from scipy import interpolate
import time
import pytest
import calendar
from datetime import datetime, timedelta

from gas_corrections_lib.src import gas_corrections

# Note: Older versions of .npy files had whole PACE image with 1710 x 1272 pixels and 197 wavelength
# i.e. reshaping arrays would look like reshape((1710, 1272, 197))
# New files only take the first 100 x 100 square of that image to save space and time

save_transmittances = False


@pytest.fixture(scope="session")
def read_PACE_data():
    l1_data = gas_corrections.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)

    return l1_data


@pytest.fixture(scope="session")
def read_MET_data():
    with Dataset('test/PACE/GMAO_MERRA2.20240411T180000.MET.nc','r') as fd:
        ozone1 = fd.variables['TO3'][:]*0.001

    with Dataset('test/PACE/GMAO_MERRA2.20240411T190000.MET.nc','r') as fd:
        ozone2 = fd.variables['TO3'][:]*0.001

    return ozone1, ozone2



def parse_timestamp():
    with h5py.File('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', 'r') as f:
        start_time = f.attrs['time_coverage_start'].decode()
        end_time = f.attrs['time_coverage_end'].decode()
    
        dt_start = time.strptime(start_time[0:19], '%Y-%m-%dT%H:%M:%S')
        dt_end = time.strptime(end_time[0:19], '%Y-%m-%dT%H:%M:%S')
        dt = time.gmtime((calendar.timegm(dt_start) + calendar.timegm(dt_end))/2)
        year = '{:04d}'.format(dt.tm_year)
        month = '{:02d}'.format(dt.tm_mon)
        day = '{:02d}'.format(dt.tm_mday)
        hour = '{:02d}'.format(dt.tm_hour)
        minute ='{:02d}'.format(dt.tm_min)
        second = '{:02d}'.format(dt.tm_sec) 

    return year, month, day, hour, minute, second


@pytest.fixture(scope="session")
def get_MET_delta_t():
    year, month, day, hour, minute, second = parse_timestamp()
    datestr = year + month + day + hour + minute + second

    l1btime = datetime.strptime(datestr,'%Y%m%d%H%M%S').timestamp()

    time1_str = datestr[0:8] + 'T' + datestr[8:10] + '0000'
    time1 = datetime.strptime(time1_str[0:8]+time1_str[9:15],'%Y%m%d%H%M%S').timestamp()
    time2 = datetime.fromtimestamp(l1btime + 3600)
    time2_str = str(time2.year) + str(time2.month).zfill(2) + str(time2.day).zfill(2) + 'T' + str(time2.hour).zfill(2) + '0000'
    ozdt = (l1btime - time1)/3600

    return ozdt


@pytest.fixture(scope="session")
def read_OCSSW_lat_lon():
    with h5py.File('test/PACE/ozone/PACE_OCI.20240411T182012.L2.ozone.nc', 'r') as f:
        OCSSW_lat = np.array(f['/navigation_data/latitude'])
        OCSSW_lon = np.array(f['/navigation_data/longitude'])

    return OCSSW_lat, OCSSW_lon


@pytest.fixture(scope="session")
def read_gas_transmittance_table():
    gas_transmittance_table = gas_corrections.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')

    return gas_transmittance_table


@pytest.fixture(scope="session")
def read_ozone_ancillary_data():
    koz = np.load('test/PACE/ozone/koz.npy')
    ozone_concentration = np.load('test/PACE/ozone/oz_concentration.npy')

    return koz, ozone_concentration


@pytest.fixture(scope="session")
def read_OCSMART_ozone_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/PACE/ozone/tg_sol_oz.npy')
    tg_sen_ocsmart = np.load('test/PACE/ozone/tg_sen_oz.npy')
    sensor_wavelengths = np.load('test/PACE/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture(scope="session")
def read_OCSSW_ozone_transmittance_benchmark_data():
    with h5py.File('test/PACE/ozone/PACE_OCI.20240411T182012.L2.ozone.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_no2_ancillary_data():
    kno2 = np.load('test/PACE/no2/k_no2.npy')
    no2_frac = np.load('test/PACE/no2/no2_frac.npy')
    no2_tropo = np.load('test/PACE/no2/no2_tropo.npy')
    no2_strat = np.load('test/PACE/no2/no2_strat.npy')

    return kno2, no2_frac, no2_tropo, no2_strat


@pytest.fixture(scope="session")
def read_OCSMART_no2_transmittance_benchmark_data():
    tg_sol_ocsmart = np.load('test/PACE/no2/tg_sol_no2.npy')
    tg_sen_ocsmart = np.load('test/PACE/no2/tg_sen_no2.npy')
    tg_sen_ocsmart = np.load('test/PACE/no2/tg_sen_no2.npy')
    sensor_wavelengths = np.load('test/PACE/OCSMART_wavelengths.npy')

    return tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths


@pytest.fixture(scope="session")
def read_OCSSW_no2_transmittance_benchmark_data():
    with h5py.File('test/PACE/no2/PACE_OCI.20240411T182012.L2.no2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_co2_transmittance_benchmark_data():
    with h5py.File('test/PACE/co2/PACE_OCI.20240411T182012.L2.co2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_co_transmittance_benchmark_data():
    with h5py.File('test/PACE/co/PACE_OCI.20240411T182012.L2.co.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_ch4_transmittance_benchmark_data():
    with h5py.File('test/PACE/ch4/PACE_OCI.20240411T182012.L2.ch4.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_n2o_transmittance_benchmark_data():
    with h5py.File('test/PACE/n2o/PACE_OCI.20240411T182012.L2.n2o.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_o2_opt2_transmittance_benchmark_data():
    with h5py.File('test/PACE/o2/PACE_OCI.20240411T182012.L2.o2.oxaband_opt_2.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
def read_OCSSW_o2_opt3_transmittance_benchmark_data():
    with h5py.File('test/PACE/o2/PACE_OCI.20240411T182012.L2.o2.oxaband_opt_3.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

    return tg_sen_ocssw, tg_sol_ocssw, wavelength_3d


@pytest.fixture(scope="session")
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
def test_ozone_OCSSW(read_OCSSW_lat_lon,
                     read_MET_data,
                     get_MET_delta_t,
                     read_ozone_ancillary_data, 
                     read_OCSMART_ozone_transmittance_benchmark_data, 
                     read_OCSSW_ozone_transmittance_benchmark_data):

    l1b_lat, l1b_lon = read_OCSSW_lat_lon
    ozone1, ozone2 = read_MET_data
    ozdt = get_MET_delta_t
    ozone_absorption_cross_section, _ = read_ozone_ancillary_data

    oz_lat = np.arange(-90.5,91.0,0.5)
    oz_lon = np.arange(-180.3125,180.625,0.625)
    oz_nline=len(oz_lat)
    oz_npixl=len(oz_lon)
    ozone=(ozone1*(1 - ozdt) + ozone2*ozdt) # interpolate in time and convert unit   
    ozmap=np.zeros((oz_nline,oz_npixl), dtype='float64')
    ozmap[1:oz_nline-1,1:oz_npixl-1]=ozone
    ozmap[:,0]=ozmap[:,oz_npixl-2]
    ozmap[:,oz_npixl-1]=ozmap[:,1]
    ozmap[0,:]=ozmap[1,:]
    ozmap[oz_nline-1,:]=ozmap[oz_nline-2,:]

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    l1_data.num_wavelengths = len(ozone_absorption_cross_section)

    # Interpolate ozone map to the L1B grid
    func = interpolate.RegularGridInterpolator((np.flip(oz_lat), oz_lon), np.flip(ozmap, 0))
    ozone_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

    ancillary_data = gas_corrections.Ancillary_Data()
    ancillary_data.ozone_absorption_cross_section = ozone_absorption_cross_section
    ancillary_data.ozone_concentration = ozone_concentration

    gas_transmittances = gas_corrections.ozone_transmittance(l1_data=l1_data, ancillary_data=ancillary_data)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_ozone_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((100, 100, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((100, 100, 197))

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
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/ozone/transmittance_comparison.png')


# @pytest.mark.skip()
def test_co2_OCSSW(read_OCSSW_co2_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.co2_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co2_transmittance_benchmark_data

    sensor_wavelengths = gas_correction_manager.l1_data.wavelengths

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
def test_co_OCSSW(read_OCSSW_co_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.co_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_co_transmittance_benchmark_data

    sensor_wavelengths = gas_correction_manager.l1_data.wavelengths

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
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/co/transmittance_comparison.png')


# @pytest.mark.skip()
def test_ch4_OCSSW(read_OCSSW_ch4_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.ch4_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_ch4_transmittance_benchmark_data

    sensor_wavelengths = gas_transmittance_manager.l1_data.wavelengths

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
def test_n2o_OCSSW(read_OCSSW_n2o_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.n2o_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_n2o_transmittance_benchmark_data

    sensor_wavelengths = gas_transmittance_manager.l1_data.wavelengths

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
    plt.ylim([0.95, 1.05])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/n2o/transmittance_comparison.png')


# @pytest.mark.skip()
def test_no2_OCSSW(read_OCSSW_lat_lon,
                   read_no2_ancillary_data, 
                   read_OCSMART_no2_transmittance_benchmark_data, 
                   read_OCSSW_no2_transmittance_benchmark_data):

    l1b_lat, l1b_lon = read_OCSSW_lat_lon
    no2_absorption_cross_section, _, _, _ = read_no2_ancillary_data
    year, month, day, hour, minute, second = parse_timestamp()

    months=range(1,13)
    nmonths=12
    #set latitude and longitude grid
    no2_frac_lat=np.arange(91,-93,-2)
    no2_frac_lon=np.arange(-181,183,2)
    no2_lat=np.arange(90.125,-90.375,-0.25)
    no2_lon=np.arange(-180.125,180.375,0.25)
    no2_frac_nline=len(no2_frac_lat)
    no2_frac_npixl=len(no2_frac_lon)
    no2_nline=len(no2_lat)
    no2_npixl=len(no2_lon)
    
    no2_total = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
    no2_tropo = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
    no2_strat = np.zeros((nmonths,no2_nline,no2_npixl), dtype='float64')
    no2_frac = np.zeros((no2_frac_nline,no2_frac_npixl), dtype='float64')
    
    no2_fname = 'test/PACE/no2/no2_climatology_v2013.hdf'
    no2_frac_fname = 'test/PACE/no2/trop_f_no2_200m.hdf'
    
    #read no2 fraction data
    f=SD(no2_frac_fname, SDC.READ)
    no2_frac[1:no2_frac_nline-1,1:no2_frac_npixl-1] = f.select('f_no2_200m')[:,:]
    no2_frac[:,0] = no2_frac[:,no2_frac_npixl-2]
    no2_frac[:,no2_frac_npixl-1] = no2_frac[:,1]
    no2_frac[0,:] = no2_frac[1,:]
    no2_frac[no2_frac_nline-1,:] = no2_frac[no2_frac_nline-2,:]
    
    # read total and tropospheric no2 data
    f=SD(no2_fname, SDC.READ)
    for i, m in enumerate(months):               
        no2_tropo[i,1:no2_nline-1,1:no2_npixl-1] = f.select('trop_no2_{:02d}'.format(m))[:,:]
        no2_tropo[i,:,0] = no2_tropo[i,:,no2_npixl-2]
        no2_tropo[i,:,no2_npixl-1] = no2_tropo[i,:,1]
        no2_tropo[i,0,:] = no2_tropo[i,1,:]
        no2_tropo[i,no2_nline-1,:] = no2_tropo[i,no2_nline-2,:]
        no2_total[i,1:no2_nline-1,1:no2_npixl-1] = f.select('tot_no2_{:02d}'.format(m))[:,:]
        no2_total[i,:,0] = no2_total[i,:,no2_npixl-2]
        no2_total[i,:,no2_npixl-1] = no2_total[i,:,1]
        no2_total[i,0,:] = no2_total[i,1,:]
        no2_total[i,no2_nline-1,:] = no2_total[i,no2_nline-2,:]
        no2_strat[i,:,:] = no2_total[i,:,:] - no2_tropo[i,:,:]
    no2_strat[no2_strat<0.0] = 0.0
    no2_total = no2_total * 1.0e15 
    no2_tropo = no2_tropo * 1.0e15
    no2_strat = no2_strat * 1.0e15

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    l1_data.num_wavelengths = len(no2_absorption_cross_section)
    
    # Interpolate no2 map to the L1B grid
    func = interpolate.RegularGridInterpolator((np.flip(no2_frac_lat), no2_frac_lon), np.flip(no2_frac, 0))
    fraction_tropospheric_no2_above_200m = func(np.array([l1b_lat, l1b_lon]).transpose())

    no2_strat = no2_strat[int(month)-1, :, :]
    func = interpolate.RegularGridInterpolator((np.flip(no2_lat), no2_lon), np.flip(no2_strat, 0))
    stratospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

    no2_tropo = no2_tropo[int(month)-1, :, :]
    func = interpolate.RegularGridInterpolator((np.flip(no2_lat), no2_lon), np.flip(no2_tropo, 0))
    tropospheric_no2_concentration = func(np.array([l1b_lat, l1b_lon]).transpose())

    ancillary_data = gas_corrections.Ancillary_Data()
    ancillary_data.no2_absorption_cross_section = no2_absorption_cross_section
    ancillary_data.fraction_tropospheric_no2_above_200m = fraction_tropospheric_no2_above_200m
    ancillary_data.tropospheric_no2_concentration = tropospheric_no2_concentration
    ancillary_data.stratospheric_no2_concentration = stratospheric_no2_concentration

    gas_transmittances = gas_corrections.no2_transmittance(l1_data=l1_data, ancillary_data=ancillary_data)

    tg_sol_ocsmart, tg_sen_ocsmart, sensor_wavelengths = read_OCSMART_no2_transmittance_benchmark_data

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape((100, 100, 197))
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape((100, 100, 197))

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
def test_o2_OCSSW_transmittance_table_option(read_OCSSW_o2_opt2_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, 
                                                          gas_transmittance_table=gas_transmittance_table, 
                                                          oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option().TRANSMITTANCE_TABLE)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_o2_opt2_transmittance_benchmark_data
    
    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'o2_table_option')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('o2_table_option')

    assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/o2_table_option/transmittance_comparison_oxaband_opt_2.png')

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([750, 800])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/o2_table_option/transmittance_comparison_oxaband_opt_2_A_band_zoom.png')


# @pytest.mark.skip()
def test_o2_OCSSW_surrounding_window_bands_option(read_OCSSW_o2_opt3_transmittance_benchmark_data):

    gas_transmittance_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_transmittance_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_transmittance_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.o2_transmittance(l1_data=l1_data, 
                                                          gas_transmittance_table=gas_transmittance_table, 
                                                          oxygen_A_band_option=gas_corrections.Oxygen_A_Band_Option().SURROUNDING_WINDOW_BANDS)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_transmittance_manager.l1_data.reflectance.shape)

    tg_sen_ocssw, tg_sol_ocssw, wavelength_3d = read_OCSSW_o2_opt3_transmittance_benchmark_data
    
    if save_transmittances:
        save_gas_transmittances(gas_transmittances, 'o2_window_bands_option')

    solar_zenith_saved, sensor_zenith_saved, total_saved = load_gas_transmittances('o2_window_bands_option')

    # assert(np.allclose(solar_zenith_saved, gas_transmittances.solar_zenith))
    # assert(np.allclose(sensor_zenith_saved, gas_transmittances.sensor_zenith))
    # assert(np.allclose(total_saved, gas_transmittances.total))

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/o2_window_bands_option/transmittance_comparison_oxaband_opt_3.png')

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :], '-r')
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :], color='orange')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_transmittance_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([750, 800])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/o2_window_bands_option/transmittance_comparison_oxaband_opt_3_A_band_zoom.png')


# @pytest.mark.skip()
def test_h2o_OCSSW(read_OCSSW_h2o_transmittance_benchmark_data):

    gas_correction_manager = gas_corrections.Gas_Correction_Manager()
    l1_data = gas_correction_manager.read_PACE_data('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', start_line=0, end_line=100, start_pixel=0, end_pixel=100)
    gas_transmittance_table = gas_correction_manager.read_gas_transmittance_table('test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc')
    gas_transmittances = gas_corrections.h2o_transmittance(l1_data=l1_data, gas_transmittance_table=gas_transmittance_table)

    tg_sen_gas_correction_lib = gas_transmittances.sensor_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)
    tg_sol_gas_correction_lib = gas_transmittances.solar_zenith.reshape(gas_correction_manager.l1_data.reflectance.shape)

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
    plt.plot(gas_correction_manager.l1_data.wavelengths[:], tg_sen_gas_correction_lib[0, 0, :], '--b')
    plt.plot(gas_correction_manager.l1_data.wavelengths[:], tg_sol_gas_correction_lib[0, 0, :], '--k')
    plt.xlim([500, 900])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transmittance')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith'])
    plt.savefig('test/PACE/h2o/transmittance_comparison.png')