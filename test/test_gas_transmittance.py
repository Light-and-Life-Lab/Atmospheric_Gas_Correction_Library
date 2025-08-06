import numpy as np

import bin.gas_transmittance as gas_transmittance

def test_ozone_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.ozone_absorption_cross_section = np.load('test/ozone_data/koz.npy')
    l1_rec.ozone_concentration = np.load('test/ozone_data/oz_concentration.npy')
    l1_rec.cos_solar_zenith = np.load('test/ozone_data/oz_csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/ozone_data/oz_csenz.npy')
    do_amf_correction = False

    t_rec = gas_transmittance.ozone_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/ozone_data/tg_sol_oz.npy')
    gas_transmittance_sensor_zenith_benchmark = np.load('test/ozone_data/tg_sen_oz.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)


def test_co_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.co_transmittance = np.load('test/co_data/t_co.npy')
    l1_rec.air_mass_factor_mixed_gases = np.load('test/co_data/amf_mixed.npy')
    l1_rec.num_amf_grid_points = np.load('test/co_data/num_airmass.npy')
    l1_rec.cos_solar_zenith = np.load('test/co_data/csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/co_data/csenz.npy')
    do_amf_correction = True

    t_rec = gas_transmittance.co_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/co_data/tg_sol_co.npy')
    gas_transmittance_total_benchmark = np.load('test/co_data/tg_co.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_total, gas_transmittance_total_benchmark)


def test_no2_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.no2_absorption_cross_section = np.load('test/no2_data/k_no2.npy')
    l1_rec.fraction_tropospheric_no2_above_200m = np.load('test/no2_data/no2_frac.npy')
    l1_rec.tropospheric_no2_concentration = np.load('test/no2_data/no2_tropo.npy')
    l1_rec.stratospheric_no2_concentration = np.load('test/no2_data/no2_strat.npy')
    l1_rec.cos_solar_zenith = np.load('test/no2_data/no2_csolz.npy')
    l1_rec.cos_sensor_zenith = np.load('test/no2_data/no2_csenz.npy')
    do_amf_correction = False

    t_rec = gas_transmittance.no2_transmittance(l1_rec, do_amf_correction)
    gas_transmittance_solar_zenith_benchmark = np.load('test/no2_data/tg_sol_no2.npy')
    gas_transmittance_sensor_zenith_benchmark = np.load('test/no2_data/tg_sen_no2.npy')

    np.testing.assert_allclose(t_rec.gas_transmittance_solar_zenith, gas_transmittance_solar_zenith_benchmark)
    np.testing.assert_allclose(t_rec.gas_transmittance_sensor_zenith, gas_transmittance_sensor_zenith_benchmark)