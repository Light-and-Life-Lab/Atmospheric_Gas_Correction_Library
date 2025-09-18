import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import interpolate
import time

import bin.gas_transmittance as gas_transmittance

# def test_ozone_transmittance():
#     ancillary_data = gas_transmittance.Ancillary_Data()
#     ancillary_data.ozone_absorption_cross_section = np.load('test/npy/ozone/koz.npy')
#     ancillary_data.ozone_concentration = np.load('test/npy/ozone/oz_concentration.npy')

#     l1_rec = gas_transmittance.L1_Record()
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(ancillary_data.ozone_absorption_cross_section)
#     do_amf_correction = False

#     t_rec = gas_transmittance.ozone_transmittance(l1_rec, ancillary_data, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/ozone/tg_sol_oz.npy')
#     gas_transmittance_sensor_zenith_benchmark = np.load('test/npy/ozone/tg_sen_oz.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)


def test_ozone_OCSSW():
    ancillary_data = gas_transmittance.Ancillary_Data()
    ancillary_data.ozone_absorption_cross_section = np.load('test/PACE/ozone/koz.npy')
    ancillary_data.ozone_concentration = np.load('test/PACE/ozone/oz_concentration.npy')

    l1_rec = gas_transmittance.L1_Record()
    l1_rec.cos_solar_zenith = np.load('test/PACE/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/PACE/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(ancillary_data.ozone_absorption_cross_section)
    do_amf_correction = False

    t_rec = gas_transmittance.ozone_transmittance(l1_rec, ancillary_data, do_amf_correction)

    shape = (1710, 1272)
    gas_transmittance_solar_zenith_benchmark = np.load('test/PACE/ozone/tg_sol_oz.npy')
    gas_transmittance_sensor_zenith_benchmark = np.load('test/PACE/ozone/tg_sen_oz.npy')

    tg_sen_ocsmart = np.flipud(gas_transmittance_solar_zenith_benchmark.reshape((1710, 1272, 197)))
    tg_sol_ocsmart = np.flipud(gas_transmittance_sensor_zenith_benchmark.reshape((1710, 1272, 197)))

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)

    with h5py.File('test/PACE/ozone/PACE_OCI.20240411T182012.L2.ozone.nc', 'r') as f:
        tg_sen_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sen'])
        tg_sol_ocssw = 5e-5*np.array(f['/geophysical_data/tg_sol'])
        wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])
        OCSSW_lat = np.array(f['/navigation_data/latitude'])
        OCSSW_lon = np.array(f['/navigation_data/longitude'])

    OCSMART_lat = np.flipud(np.load('test/PACE/l1b_lat.npy').reshape((1710, 1272)))
    OCSMART_lon = np.flipud(np.load('test/PACE/l1b_lon.npy').reshape((1710, 1272)))

    tg_sen_gas_correction_lib = np.flipud(t_rec.gas_transmittance_sensor_zenith.reshape((1710, 1272, 197)))
    tg_sol_gas_correction_lib = np.flipud(t_rec.gas_transmittance_solar_zenith.reshape((1710, 1272, 197)))

    sensor_wavelengths = np.load('test/PACE/OCSMART_wavelengths.npy')

    # OCSMART_first_pixel_idx = 2173848 # Index of the location (6.67, -93.1), which is the first pixel in the OCSSW grid, but not in the OCSMART grid

    plt.figure()
    plt.plot(wavelength_3d, tg_sen_ocssw[0, 0, :])
    plt.plot(wavelength_3d, tg_sol_ocssw[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_gas_correction_lib[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sol_gas_correction_lib[0, 0, :])
    plt.plot(sensor_wavelengths, tg_sen_ocsmart[0, 0, :], '--b')
    plt.plot(sensor_wavelengths, tg_sol_ocsmart[0, 0, :], '--k')
    plt.legend(['OCSSW Sensor Zenith', 'OCSSW Solar Zenith', 'Gas Lib Sensor Zenith', 'Gas Lib Solar Zenith', 'OCSMART Sensor Zenith', 'OCSMART Solar Zenith'])
    plt.savefig('test/PACE/ozone/transmittance_comparison.png')


# def test_co2_transmittance():
#     amf_table = gas_transmittance.Air_Mass_Factor_Lookup_Table()
#     amf_table.co2_transmittance = np.load('test/npy/co2/t_co2.npy')
#     amf_table.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
#     amf_table.num_amf_grid_points = np.load('test/npy/num_airmass.npy')

#     l1_rec = gas_transmittance.L1_Record()
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(amf_table.co2_transmittance[:, 0])
#     do_amf_correction = True

#     t_rec = gas_transmittance.co2_transmittance(l1_rec, amf_table, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/co2/tg_sol_co2.npy')
#     gas_transmittance_total_benchmark = np.load('test/npy/co2/tg_co2.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


# def test_co_transmittance():
#     amf_table = gas_transmittance.Air_Mass_Factor_Lookup_Table()
#     amf_table.co_transmittance = np.load('test/npy/co/t_co.npy')
#     amf_table.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
#     amf_table.num_amf_grid_points = np.load('test/npy/num_airmass.npy')

#     l1_rec = gas_transmittance.L1_Record()
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(amf_table.co_transmittance[:, 0])
#     do_amf_correction = True

#     t_rec = gas_transmittance.co_transmittance(l1_rec, amf_table, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/co/tg_sol_co.npy')
#     gas_transmittance_total_benchmark = np.load('test/npy/co/tg_co.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


# def test_ch4_transmittance():
#     amf_table = gas_transmittance.Air_Mass_Factor_Lookup_Table()
#     amf_table.ch4_transmittance = np.load('test/npy/ch4/t_ch4.npy')
#     amf_table.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
#     amf_table.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    
#     l1_rec = gas_transmittance.L1_Record()
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(amf_table.ch4_transmittance[:, 0])
#     do_amf_correction = True

#     t_rec = gas_transmittance.ch4_transmittance(l1_rec, amf_table, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/ch4/tg_sol_ch4.npy')
#     gas_transmittance_total_benchmark = np.load('test/npy/ch4/tg_ch4.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


# def test_n2o_transmittance():
#     amf_table = gas_transmittance.Air_Mass_Factor_Lookup_Table()
#     amf_table.n2o_transmittance = np.load('test/npy/n2o/t_n2o.npy')
#     amf_table.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
#     amf_table.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    
#     l1_rec = gas_transmittance.L1_Record()
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(amf_table.n2o_transmittance[:, 0])
#     do_amf_correction = True

#     t_rec = gas_transmittance.n2o_transmittance(l1_rec, amf_table, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/n2o/tg_sol_n2o.npy')
#     gas_transmittance_total_benchmark = np.load('test/npy/n2o/tg_n2o.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


# def test_no2_transmittance():
#     ancillary_data = gas_transmittance.Ancillary_Data()
#     ancillary_data.no2_absorption_cross_section = np.load('test/npy/no2/no2_absorption_cross_section.npy')
#     ancillary_data.fraction_tropospheric_no2_above_200m = np.load('test/npy/no2/no2_frac.npy')
#     ancillary_data.tropospheric_no2_concentration = np.load('test/npy/no2/no2_tropo.npy')
#     ancillary_data.stratospheric_no2_concentration = np.load('test/npy/no2/no2_strat.npy')
    
#     l1_rec = gas_transmittance.L1_Record() 
#     l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
#     l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
#     l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
#     l1_rec.num_wavelengths = len(ancillary_data.no2_absorption_cross_section)
#     do_amf_correction = False

#     t_rec = gas_transmittance.no2_transmittance(l1_rec, ancillary_data, do_amf_correction)
#     gas_transmittance_solar_zenith_benchmark = np.load('test/npy/no2/tg_sol_no2.npy')
#     gas_transmittance_sensor_zenith_benchmark = np.load('test/npy/no2/tg_sen_no2.npy')

#     np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)


# def test_o2_transmittance():
#     amf_table = gas_transmittance.Air_Mass_Factor_Lookup_Table()

#     with h5py.File("test/PACE/oci_gas_transmittance_cia_amf_v3.2.nc", 'r') as f:
#         o2_transmittance = np.array(f['oxygen_transmittance'])
#         amf_table.air_mass_factor_mixed_gases = np.array(f['air_mass_factor_mixed'])
#         amf_table.num_amf_grid_points = np.array(f['n_air_mass_factor']).size
#         gas_transmittance_wavelengths = np.array(f['wavelength'])

#     l1_rec = gas_transmittance.L1_Record()
#     start_line = 0
#     start_pixel = 0
#     end_line = 100
#     end_pixel = 100

#     with h5py.File('test/PACE/PACE_OCI.20240411T182012.L1B.V3.nc', 'r') as f:
#         l1_rec.cos_solar_zenith = np.cos(np.deg2rad(np.array(f['/geolocation_data/solar_zenith'][start_line:end_line, start_pixel:end_pixel])))
#         l1_rec.cos_sensor_zenith = np.cos(np.deg2rad(np.array(f['/geolocation_data/sensor_zenith'][start_line:end_line, start_pixel:end_pixel])))
#         blue_wavelengths = np.array(f['/sensor_band_parameters/blue_wavelength'][1:])
#         red_wavelengths = np.array(f['/sensor_band_parameters/red_wavelength'][3:])
#         sensor_wavelengths = np.zeros(len(blue_wavelengths) + len(red_wavelengths))
#         sensor_wavelengths[0:len(blue_wavelengths)] = blue_wavelengths
#         sensor_wavelengths[len(blue_wavelengths):] = red_wavelengths

#         blue_rhot = np.array(f['/observation_data/rhot_blue'][1:, start_line:end_line, start_pixel:end_pixel])
#         red_rhot = np.array(f['/observation_data/rhot_red'][3:, start_line:end_line, start_pixel:end_pixel])
#         assert(blue_rhot.shape[1] == red_rhot.shape[1])
#         assert(blue_rhot.shape[2] == red_rhot.shape[2])

#         rhot = np.zeros((blue_rhot.shape[0] + red_rhot.shape[0], blue_rhot.shape[1], blue_rhot.shape[2]))
#         rhot[0:len(blue_wavelengths), :, :] = blue_rhot
#         rhot[len(blue_wavelengths):, :, :] = red_rhot
#         rhot = np.rollaxis(rhot, 0, 3)
#         reflectance = rhot/np.pi*l1_rec.cos_sensor_zenith[:, :, None]

#         assert(len(sensor_wavelengths) == reflectance.shape[2])

#         l1_rec.wavelengths = sensor_wavelengths
#         l1_rec.num_pixels = len(l1_rec.cos_solar_zenith.flatten())
#         l1_rec.num_wavelengths = len(l1_rec.wavelengths)
#         do_amf_correction = False

#     with h5py.File('test/PACE/PACE_OCI.20240411T182012.L2.nc', 'r') as f:
#         tg_sen_benchmark = 5e-5*np.array(f['/geophysical_data/tg_sen'])
#         tg_sol_benchmark = 5e-5*np.array(f['/geophysical_data/tg_sol'])
#         wavelength_3d = np.array(f['/sensor_band_parameters/wavelength_3d'])

#     plt.figure()
#     plt.plot(wavelength_3d, tg_sen_benchmark[0, 0, :])
#     plt.plot(wavelength_3d, tg_sol_benchmark[0, 0, :])
#     plt.savefig('test/PACE/benchmark_transmittances.png')

#     F0_gt_wavelengths = [112.329, 92.673, 85.208, 82.100, 80.692, 86.329, 95.925, 101.478, 101.788, 98.128, 93.719, 95.656, 97.642, 95.100, 93.133, 95.137, 99.255, 103.093, 99.177, 93.290, 97.047, 105.519, 114.651, 119.597, 116.077, 107.558, 109.618, 119.291, 112.055, 97.243, 97.233, 107.436, 112.781, 106.594, 105.608, 126.993, 156.889, 170.586, 169.248, 169.343, 173.391, 177.938, 178.901, 176.584, 175.633, 174.629, 171.281, 162.301, 153.935, 161.212, 174.113, 178.181, 182.051, 190.560, 195.324, 199.413, 205.297, 205.224, 204.442, 206.345, 207.477, 208.190, 206.088, 203.623, 203.584, 205.278, 207.110, 208.286, 209.299, 206.511, 195.706, 190.027, 195.623, 199.210, 200.467, 199.674, 195.068, 193.092, 195.670, 197.349, 196.529, 193.713, 186.203, 179.120, 181.903, 188.641, 188.972, 188.225, 191.667, 192.903, 192.064, 190.857, 188.243, 188.456, 189.880, 189.660, 189.365, 189.610, 188.825, 185.808, 184.117, 184.857, 184.926, 184.325, 184.315, 185.745, 185.854, 184.273, 184.355, 184.525, 182.317, 178.484, 177.707, 179.744, 179.969, 178.360, 176.334, 176.103, 176.462, 174.718, 172.349, 170.259, 167.958, 167.776, 168.949, 167.509, 165.836, 166.449, 165.916, 164.217, 163.716, 163.528, 162.046, 161.215, 160.637, 160.358, 160.231, 159.892, 159.257, 158.725, 158.332, 158.267, 157.786, 155.176, 151.696, 148.458, 147.693, 149.665, 152.508, 154.888, 155.421, 155.288, 154.996, 154.702, 154.397, 154.032, 153.450, 152.850, 152.293, 151.891, 151.637, 151.432, 151.158, 150.793, 150.431, 149.908, 149.250, 148.497, 147.875, 147.506, 147.418, 147.395, 147.151, 146.716, 146.142, 145.677, 145.179, 144.677, 144.075, 143.268, 142.395, 141.534, 141.038, 140.862, 140.959, 140.985, 140.752, 140.295, 139.692, 139.124, 138.487, 137.908, 137.372, 136.140, 134.717, 134.219, 134.304, 133.385, 132.080, 131.723, 131.437, 130.196, 128.322, 127.705, 127.371, 127.484, 127.690, 127.841, 127.726, 127.292, 126.747, 126.200, 125.866, 125.688, 125.509, 125.165, 124.782, 124.408, 124.097, 123.926, 123.728, 123.220, 122.490, 121.679, 121.056, 120.707, 120.507, 120.407, 120.035, 119.765, 119.542, 119.096, 118.633, 118.122, 117.631, 117.233, 116.445, 114.587, 113.353, 113.677, 113.456, 112.543, 111.714, 110.849, 110.234, 110.170, 109.946, 108.613, 106.799, 106.192, 106.574, 106.537, 105.791, 104.368, 103.511, 103.470, 103.245, 102.549, 101.870, 100.769, 98.572, 94.626, 92.347, 95.812, 98.548, 97.907, 94.326, 91.864, 94.043, 95.080, 94.395, 94.003, 93.393, 93.138, 92.918, 92.343, 92.274, 91.902, 91.528, 81.981, 67.021, 44.511, 35.562, 23.500, 9.111, 7.397] # TODO: Need to figure out how to get real values for F0 in OC-SMART. OCSSW reads them in from sensor_info.dat files
#     f = interpolate.interp1d(gas_transmittance_wavelengths, F0_gt_wavelengths)
#     F0_sensor_wavelengths = f(sensor_wavelengths)

#     f = interpolate.interp1d(gas_transmittance_wavelengths, o2_transmittance, axis = 0)
#     o2_transmittance_sensor_wavelengths = f(sensor_wavelengths)

#     l1_rec.Lt = reflectance
#     l1_rec.F0 = F0_sensor_wavelengths
#     amf_table.o2_transmittance = o2_transmittance_sensor_wavelengths

#     oxygen_A_band_option = gas_transmittance.Oxygen_A_Band_Option.NO_AMF_CORRECTION

#     start_time = time.perf_counter()
#     t_rec = gas_transmittance.o2_transmittance(l1_rec, amf_table, do_amf_correction, oxygen_A_band_option)
#     end_time = time.perf_counter()

#     plt.plot(sensor_wavelengths, t_rec.gas_transmittance_sensor_zenith[0, :])
#     plt.plot(sensor_wavelengths, t_rec.gas_transmittance_solar_zenith[0, :])
#     plt.savefig('test/PACE/benchmark_transmittances.png')

#     cpp_time_elapsed = end_time - start_time
#     print("C++ O2 time: ", cpp_time_elapsed)


#     # gas_transmittance_solar_zenith_benchmark = np.load('test/npy/n2o/tg_sol_n2o.npy')
#     # gas_transmittance_total_benchmark = np.load('test/npy/n2o/tg_n2o.npy')

#     # np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
#     # np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)



#     # TODO: Implement this test
#     def test_h2o_transmittance():
#         pass