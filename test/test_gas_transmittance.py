import numpy as np

import bin.gas_transmittance as gas_transmittance

def test_ozone_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.ozone_absorption_cross_section = np.load('test/npy/ozone/koz.npy')
    l1_rec.ozone_concentration = np.load('test/npy/ozone/oz_concentration.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.ozone_absorption_cross_section)
    do_amf_correction = False

    t_rec = gas_transmittance.ozone_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/ozone/tg_sol_oz.npy')
    gas_transmittance_sensor_zenith_benchmark = np.load('test/npy/ozone/tg_sen_oz.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)


def test_co2_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.co2_transmittance = np.load('test/npy/co2/t_co2.npy')
    l1_rec.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
    l1_rec.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.co2_transmittance[:, 0])
    do_amf_correction = True

    t_rec = gas_transmittance.co2_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/co2/tg_sol_co2.npy')
    gas_transmittance_total_benchmark = np.load('test/npy/co2/tg_co2.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


def test_co_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.co_transmittance = np.load('test/npy/co/t_co.npy')
    l1_rec.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
    l1_rec.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.co_transmittance[:, 0])
    do_amf_correction = True

    t_rec = gas_transmittance.co_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/co/tg_sol_co.npy')
    gas_transmittance_total_benchmark = np.load('test/npy/co/tg_co.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


def test_ch4_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.ch4_transmittance = np.load('test/npy/ch4/t_ch4.npy')
    l1_rec.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
    l1_rec.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.ch4_transmittance[:, 0])
    do_amf_correction = True

    t_rec = gas_transmittance.ch4_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/ch4/tg_sol_ch4.npy')
    gas_transmittance_total_benchmark = np.load('test/npy/ch4/tg_ch4.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


def test_n2o_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.n2o_transmittance = np.load('test/npy/n2o/t_n2o.npy')
    l1_rec.air_mass_factor_mixed_gases = np.load('test/npy/amf_mixed.npy')
    l1_rec.num_amf_grid_points = np.load('test/npy/num_airmass.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.n2o_transmittance[:, 0])
    do_amf_correction = True

    t_rec = gas_transmittance.n2o_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/n2o/tg_sol_n2o.npy')
    gas_transmittance_total_benchmark = np.load('test/npy/n2o/tg_n2o.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


def test_no2_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.no2_absorption_cross_section = np.load('test/npy/no2/k_no2.npy')
    l1_rec.fraction_tropospheric_no2_above_200m = np.load('test/npy/no2/no2_frac.npy')
    l1_rec.tropospheric_no2_concentration = np.load('test/npy/no2/no2_tropo.npy')
    l1_rec.stratospheric_no2_concentration = np.load('test/npy/no2/no2_strat.npy')
    l1_rec.cos_solar_zenith = np.load('test/npy/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/npy/csenz.npy')
    l1_rec.num_pixels = len(l1_rec.cos_solar_zenith)
    l1_rec.num_wavelengths = len(l1_rec.no2_absorption_cross_section)
    do_amf_correction = False

    t_rec = gas_transmittance.no2_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/npy/no2/tg_sol_no2.npy')
    gas_transmittance_sensor_zenith_benchmark = np.load('test/npy/no2/tg_sen_no2.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)