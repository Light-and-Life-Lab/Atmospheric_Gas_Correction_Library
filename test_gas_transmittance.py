import numpy as np
import cppimport
funcs = cppimport.imp("gas_transmittance")

import gas_transmittance

if __name__ == "__main__":
    k_oz = np.array([1, 2, 3])
    l1b_oz = np.ones(4)
    l1b_solz = np.array([90, 60, 30, 0])
    l1b_senz = np.array([90, 60, 30, 0])
    do_amf_correction = False

    tg_sol, tg_sen, tg = gas_transmittance.ozone_transmittance(k_oz, l1b_oz, l1b_solz, l1b_senz, do_amf_correction)
    print(tg_sol)
    print(tg_sen)
    print(tg)