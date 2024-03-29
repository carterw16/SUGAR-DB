import numpy as np
import pickle

def extract_results(multi) -> dict:
    """ Extracts the following grid state variables from the Multi object for
    each period in the final converged epoch:
        - node: electrical bus - V = voltage (complex phasor in volts)
        - load: constant power PQ loads (kW) - includes generation of PV and wind
          as negative PQ loads
        - ohline: overhead lines - S_oh = power flow (kVA)
        - ugline: underground liens - S_ug = power flow (kVA) 
        - slack: slack generator - P_g = generation power (kW)
        - battery has 3 phases, following variables in length 3 list [PhaseA, PhaseB, PhaseC]: 
            B = state of charge (% of total)  
            P_ch = charge power, 
            P_d = discharge power (kW)
        - total load (kW)
    Multi-period Metrics:
        - total cost C = $
        - solve time (s)
    """

    num_bats = Multi.num_batteries
    num_lines = Multi.infeasibility.casedata.ohlines 
    results = {'C_tot': 0, 'Solve_time': 0, 'B': np.zeros((Multi.periods, num_bats, 3)), 
            'P_ch': np.zeros((Multi.periods, num_bats, 3)), 
            'P_d': np.zeros((Multi.periods, num_bats, 3)), 
            'P_g': np.zeros((Multi.periods)), 'S_ug': [], 'S_oh': [], 'V': []
            }

    for period in Multi.periods:
        Vsol = Multi.Vsol_mp[period]
        for idx in num_bats:
            bat = Multi.battery_sources[idx]

            B = Vsol[bat.Bt_nodes]
            P_ch = Vsol[[bat.nodeA_p_plus,bat.nodeB_p_plus,bat.nodeC_p_plus]].ravel()
            P_d = Vsol[[bat.nodeA_p_minus,bat.nodeB_p_minus,bat.nodeC_p_minus]].ravel()

            results['B'][idx][period] = B
            results['P_ch'][idx][period] = P_ch
            results['P_d'][idx][period] = P_d

        for line in Multi.infeasibility.


    
    return results
