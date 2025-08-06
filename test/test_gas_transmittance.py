import numpy as np

import bin.gas_transmittance as gas_transmittance

def test_ozone_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.k_oz = np.load('test/ozone_data/koz.npy')
    l1_rec.l1b_oz = np.load('test/ozone_data/oz_concentration.npy')
    l1_rec.l1b_csolz = np.load('test/ozone_data/oz_csolz.npy')
    l1_rec.l1b_csenz = np.load('test/ozone_data/oz_csenz.npy')
    do_amf_correction = False

    t_rec = gas_transmittance.ozone_transmittance(l1_rec, do_amf_correction)
    tg_sol_benchmark = np.load('test/ozone_data/tg_sol_oz.npy')
    tg_sen_benchmark = np.load('test/ozone_data/tg_sen_oz.npy')

    np.testing.assert_allclose(t_rec.tg_sol, tg_sol_benchmark)
    np.testing.assert_allclose(t_rec.tg_sen, tg_sen_benchmark)


def test_co_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.t_co = np.load('test/co_data/t_co.npy')
    l1_rec.amf_mixed = np.load('test/co_data/amf_mixed.npy')
    l1_rec.num_airmass = np.load('test/co_data/num_airmass.npy')
    l1_rec.l1b_csolz = np.load('test/co_data/csolz.npy')
    l1_rec.l1b_csenz = np.load('test/co_data/csenz.npy')
    do_amf_correction = True

    t_rec = gas_transmittance.co_transmittance(l1_rec, do_amf_correction)
    tg_sol_benchmark = np.load('test/co_data/tg_sol_co.npy')
    tg_benchmark = np.load('test/co_data/tg_co.npy')

    np.testing.assert_allclose(t_rec.tg_sol, tg_sol_benchmark)
    np.testing.assert_allclose(t_rec.tg, tg_benchmark)


def test_no2_transmittance():
    l1_rec = gas_transmittance.L1_Record()
    l1_rec.k_no2 = np.load('test/no2_data/k_no2.npy')
    l1_rec.l1b_no2_frac = np.load('test/no2_data/no2_frac.npy')
    l1_rec.l1b_no2_tropo = np.load('test/no2_data/no2_tropo.npy')
    l1_rec.l1b_no2_strat = np.load('test/no2_data/no2_strat.npy')
    l1_rec.l1b_csolz = np.load('test/no2_data/no2_csolz.npy')
    l1_rec.l1b_csenz = np.load('test/no2_data/no2_csenz.npy')
    do_amf_correction = False

    t_rec = gas_transmittance.no2_transmittance(l1_rec, do_amf_correction)
    tg_sol_benchmark = np.load('test/no2_data/tg_sol_no2.npy')
    tg_sen_benchmark = np.load('test/no2_data/tg_sen_no2.npy')

    np.testing.assert_allclose(t_rec.tg_sol, tg_sol_benchmark)
    np.testing.assert_allclose(t_rec.tg_sen, tg_sen_benchmark)