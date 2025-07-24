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

    assert np.allclose(tg_sol, tg_sol_benchmark)
    assert np.allclose(tg_sen, tg_sen_benchmark)

