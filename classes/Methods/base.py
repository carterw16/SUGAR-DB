import numpy as np
import pickle

from classes.CurrentMeasurements import create_current_meas
from classes.Homotopy import Homotopy
from classes.Limiting import Limiting
from classes.OperationalLimits.VoltageUnbalance import VoltageUnbalance

from lib.save_output import save_output

class Methods():
    def __init__(self, testcase, casedata, features, settings, path_to_output, node_key, node_index_):
        self.features = features
        self.settings = settings
        self.path_to_output = path_to_output
        self.case_name = testcase

        self.NR_tolerance = settings["Tolerance"]

        # Initialization Settings
        self.warm_start = settings["Warm Start"]["initialize"]
        self.load_factor = settings["Load Factor"]
        self.max_iter = settings["Max Iters"]

        # Operational Settings
        if settings["voltage unbalance settings"]["enforce voltage unbalance"]:
            self.voltage_unbalance_eqn_enabled = True
            self.voltage_unbalance_settings = settings["voltage unbalance settings"]
            self.voltage_unbalance = VoltageUnbalance(self.voltage_unbalance_settings)
        else:
            self.voltage_unbalance_eqn_enabled = False
        
        if settings["voltage bound settings"]["enforce voltage bounds"]:
            self.voltage_bounds_enabled = True
            self.voltage_bound_settings = settings["voltage bound settings"]
            self.vmax_pu = self.voltage_bound_settings["Vmag Bound Max pu"]
            self.vmin_pu = self.voltage_bound_settings["Vmag Bound Min pu"]
        else:
            self.voltage_bounds_enabled = False
            self.vmax_pu = 1.05
            self.vmin_pu = 0.95

        # Output Settings
        self.save_output_conditions_flag = False
        self.draw_network_maps = True
        self.save_results = settings["Save File"]
        
        self.node = casedata.node
        self.load = casedata.load
        self.slack = casedata.slack
        self.triplex_load = casedata.triplex_load
        self.xfmr = casedata.xfmr
        self.line_oh = casedata.ohline
        self.line_ug = casedata.ugline
        self.line_tplx = casedata.triplex_line
        self.capacitor = casedata.capacitor
        self.regulator = casedata.regulator
        self.switch = casedata.switch
        self.fuse = casedata.fuse
        self.reactor = casedata.reactor
        self.ibdg = casedata.ibdg
        self.simulation_stats = casedata.stats
        self.voltage_bounds = casedata.voltage_bounds

        self.node_index_ = node_index_
        self.node_key = node_key

        """---------------------- Initialize Heuristic Settings -------------------------"""        
        # Homotopy
        if settings["Homotopy"]:
            self.homotopy_enabled = True
            g_homotopy = settings["G_homotopy"]
            b_homotopy = settings["B_homotopy"]
            h_factor_init = 1
            homotopy_tolerance = 1e-3
            h_step_size_init = 0.1   
            self.homotopy = Homotopy(homotopy_tolerance, g_homotopy, b_homotopy, h_factor_init, h_step_size_init) 
        else:
            self.homotopy = None
            self.homotopy_enabled = False
       
        # Voltage Limiting
        if settings["Voltage Limiting"]:
            voltage_limiting_dict = {'maxstep': 0.1, 'minstep': -0.1, 
                                    'Vr Limit': 2, 'Vi Limit': 2, 'dual': True, 'Vr Norm': None, 'Vi Norm': None}
            self.voltage_limiting = Limiting('voltage', voltage_limiting_dict)
            self.voltage_limiting_enabled = True
        else:
            self.voltage_limiting_enabled = False
        
        # Diode Limiting 
        if settings["Diode Limiting"]:
            self.diode_limiting_enabled = True
        else:
            self.diode_limiting_enabled = False

        enable_CM = self.features['Current Meas']
        self.curr_meas = []
        if enable_CM:
            self.activate_current_measures(enable_CM, self.node_index_, self.node_key)

        """---------------------- Initialize Operational Settings -------------------------"""    
        if settings["voltage unbalance settings"]["enforce voltage unbalance"]:
            self.voltage_unbalance_eqn_enabled = True
            self.voltage_unbalance_settings = settings["voltage unbalance settings"]

    
    def almostEqual(d1, d2, epsilon= 10**-7):
            return abs(d2 - d1) < epsilon 
    
    def activate_current_measures(self, enable_CM, node_index_, node_key):
        self.curr_meas = create_current_meas(enable_CM, self.line_oh, self.line_ug,
                                             self.line_tplx)

        for ele in self.curr_meas:
            node_index_ = ele.assign_nodes(node_index_, node_key, self.node)

        self.node_index_ = node_index_ 

    def save_output_conditions(self, node_name, voltage_name, Vsol):
        name_node = "Initial_Conditions/" + self.case_name + node_name + str(self.load_factor) + ".pkl"
        with open(name_node, "wb") as f:
            pickle.dump(self.node, f)
                    
        name_V =  "Initial_Conditions/" + self.case_name + voltage_name + str(self.load_factor) + ".npy"
        np.save(name_V, Vsol)
    
    def write_results(self, V, iteration_count, simulation_stats, dual_info_ABC, dual_info_tplx):

        for ele in self.node:
            ele.print_slack = True
            ele.calcMagAng(V, False)

        if self.regulator:
            for ele in self.regulator:
                ele.update_values(V)

        enable_IBDGs = True if self.ibdg else False
        if self.ibdg:
            for ele in self.ibdg:
                ele.calc_ibdg_outputs(self.node, self.node_key, V)

        for ele in self.curr_meas:
            ele.calc_currents(V)
            ele.calc_powers(V)

        for ele in self.xfmr:
            ele.calc_currents(V)

        for ele in self.fuse:
            ele.calc_currents(V)

        for ele in self.switch:
            ele.calc_currents(V)

        outputs = [
            self.load, self.ibdg, self.regulator, self.curr_meas, self.xfmr, self.fuse, self.switch,
            self.triplex_load, dual_info_ABC, dual_info_tplx
        ]

        # Write csv file with results
        save_output(self.path_to_output, self.case_name, self.node, outputs, iteration_count,
                    simulation_stats, self.settings, enable_IBDGs, self.load_factor, self.stamp_dual)

       
