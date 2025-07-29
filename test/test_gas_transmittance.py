import numpy as np

import bin.gas_transmittance as gas_transmittance

def test_ozone_transmittance():
    k_oz = np.load('test/ozone_data/koz.npy')
    l1b_oz = np.load('test/ozone_data/oz_concentration.npy')
    l1b_csolz = np.load('test/ozone_data/oz_csolz.npy')
    l1b_csenz = np.load('test/ozone_data/oz_csenz.npy')
    do_amf_correction = False

    tg_sol, tg_sen, tg = gas_transmittance.ozone_transmittance(k_oz, l1b_oz, l1b_csolz, l1b_csenz, do_amf_correction)
    tg_sol_benchmark = np.load('test/ozone_data/tg_sol_oz.npy')
    tg_sen_benchmark = np.load('test/ozone_data/tg_sen_oz.npy')

    np.testing.assert_allclose(tg_sol, tg_sol_benchmark)
    np.testing.assert_allclose(tg_sen, tg_sen_benchmark)


def test_no2_transmittance():
    k_no2 = np.load('test/no2_data/k_no2.npy')
    l1b_no2_tropo = np.load('test/no2_data/no2_tropo.npy')
    l1b_no2_strat = np.load('test/no2_data/no2_strat.npy')
    l1b_no2_frac = np.load('test/no2_data/no2_frac.npy')
    l1b_csolz = np.load('test/no2_data/no2_csolz.npy')
    l1b_csenz = np.load('test/no2_data/no2_csenz.npy')
    do_amf_correction = False

    tg_sol, tg_sen, tg = gas_transmittance.no2_transmittance(k_no2, l1b_no2_frac, l1b_no2_tropo, l1b_no2_strat, l1b_csolz, l1b_csenz, do_amf_correction)
    tg_sol_benchmark = np.load('test/no2_data/tg_sol_no2.npy')
    tg_sen_benchmark = np.load('test/no2_data/tg_sen_no2.npy')

    np.testing.assert_allclose(tg_sol, tg_sol_benchmark)
    np.testing.assert_allclose(tg_sen, tg_sen_benchmark)