
import logging
import math
import numpy as np
import time

from itertools import count
from scipy.sparse.linalg import spsolve
from termcolor import colored

from classes.Limiting import Limiting
from classes.Map import Map
from classes.Methods.base import Methods
from classes.Nodes import Nodes
from classes.OperationalLimits.operational_limits_methods import operational_limits_diode_limiting, calc_voltage_violations
from classes.OperationalLimits.VoltageLimits import VoltageLimits

from lib.check_zero_rows import check_zero_rows
from lib.dual_analysis import dual_analysis
from lib.initialize import initialize
from lib.initialize_warm_start import initialize_warm_start, initialize_from_amrit_solution
from lib.normalize_if import normalize_if
from lib.stamp_linear import stamp_linear
from lib.stamp_nonlinear import stamp_nonlinear
from lib.matrix_conversion import matrix_conversion
from lib.read_open_dss_results import read_open_dss_results

import ipdb


class InfeasibilityAnalysis(Methods):
    def __init__(self, testcase, casedata, features, settings, path_to_output, node_key, node_index_):
        Methods.__init__(self, testcase, casedata, features, settings, path_to_output, node_key, node_index_)

        # Infeasibility Settings
        self.infeasibility_sources = casedata.infeasibility_sources
        self.infeas_settings = settings['infeas settings']
        self.no_print = self.infeas_settings['no print']
        self.stamp_dual = settings["Stamp Dual"]
        self.stamp_neutral_infeas_source_flag = self.infeas_settings["neutral infeas source"]
        self.use_Vinit = False

        # Create Battery Sources
        self.battery_sources = casedata.battery_sources

        try:
            self.source_type = self.infeas_settings['source type']
        except KeyError:
            print("Using default infeasibility source type")
            self.source_type = 'current'

        try:
            self.obj = self.infeas_settings['obj']
        except KeyError:
            if self.stamp_dual:
                print("Using default objective minimization function")
                self.obj = 'L2'
            else:
                print('Defaulting to power flow')
                self.obj = 'PF'
        try:
            self.obj_type = self.infeas_settings['obj type']
        except KeyError:
            print("Using default objective function parameter")
            self.obj_type = 'power'

        try:
            self.stamp_slack_flag = self.infeas_settings["stamp slack bus"]
        except KeyError:
            print("Using default Slack bus settings for infeasibility sources")
            self.stamp_slack_flag = True

        try:
            self.report_if = self.infeas_settings['report infeasibility']
        except KeyError:
            self.report_if = False

        try:
            self.cs_eps = float(self.infeas_settings['comp slack tol'])
        except KeyError:
            print("Using default cs eps value for complementary slackness")
            self.cs_eps = 1e-6

        try:
            self.stamp_tplx_infeas_sources_flag = self.infeas_settings["triplex sources"]
        except KeyError:
            self.stamp_tplx_infeas_sources_flag = True

        # Output Settings
        try:
            self.normalize_if = self.infeas_settings['normalize mapped infeas var']
        except KeyError:
            self.normalize_if = False

        try:
            self.draw_infeasibility_maps = self.infeas_settings['draw infeasibility maps']
        except KeyError:
            self.draw_infeasibility_maps = False

    def calc_residuals(self, V_out, Ylin, Jlin, node_key):
        # Initializing the vector which houses the residual equations
        res_eqn = np.zeros_like(V_out)

        violation = 0

        # Stationarity Constraints (L2 norm)
        # dL/dif = 0 (located in if node) [all elements for this should be linear and this in Ylin, Jlin]

        # Primal Feasibility (Equality Constraints) [some elements in this equation are nonlinear, so these components
        # are directly added to res_eqn for load and triplex load]
        # eq constr = 0 (located in V nodes)
        for ele in range(len(self.load)):
            res_eqn = self.load[ele].calc_residual(V_out, node_key, self.node, res_eqn, self.load_factor)

        for ele in range(len(self.triplex_load)):
            res_eqn = self.triplex_load[ele].calc_residual(V_out, node_key, self.node, res_eqn, self.load_factor)

        if self.voltage_bounds:
            for ele in self.voltage_bounds:
                ele.calc_residuals(V_out, res_eqn)

        if self.voltage_unbalance_eqn_enabled:
            for ele in self.node:
                if ele.phases == 7 or ele.phases == 15:
                    res_eqn = self.voltage_unbalance.calc_residuals(V_out, ele, res_eqn)

        if self.source_type == 'PQ' or self.source_type == 'Q':
            for ele in self.infeasibility_sources:
                infeas_connected_node = self.node_key[ele.connected_node_name]
                res_eqn = ele.calc_residuals(self.node[infeas_connected_node], V_out, res_eqn, self.cs_eps)

        if self.source_type == 'B' or self.source_type == 'GB':
            print('Need to add residuals for G/B sources')

        if self.obj == 'L1' and self.source_type == 'current':
            prim_feas_flag = np.zeros((len(self.infeasibility_sources),1))
            dual_feas_flag = np.zeros((len(self.infeasibility_sources),1))
            for ele in range(len(self.infeasibility_sources)):
                infeas_connected_node = self.node_key[self.infeasibility_sources[ele].connected_node_name]
                res_eqn, prim_flag, dual_flag = self.infeasibility_sources[ele].calc_residuals_L1(V_out, self.node[infeas_connected_node], self.cs_eps, res_eqn)
                prim_feas_flag[ele] = prim_flag
                dual_feas_flag[ele] = dual_flag


        # calc residuals for battery
        for ele in self.battery_sources:
            bat_node = self.node_key[ele.connected_node_name]
            res_eqn = ele.calc_residuals(self.node[bat_node], V_out, res_eqn, self.cs_eps)

        # If res_eqn = 0, then stationarity, primal feasibility (equality), and complementary slackness are met
        res_eqn += (Ylin@V_out - Jlin)

        if any(abs(res_eqn[np.append(Nodes.Vr_index, Nodes.Vi_index)]) > self.NR_tolerance):
            print('Primal feasibility (equality) condition violated')
            violation = 1

        if self.stamp_dual == True:
            if any(abs(res_eqn[Nodes.L_index]) > self.NR_tolerance):
                print('Stationarity condition (dL/dVr) violated')
                violation = 1

            if self.source_type == 'current':
                if self.obj == "L1":
                    if np.any(abs(res_eqn[np.append(np.append(Nodes.if_r_plus_index, Nodes.if_i_plus_index), np.append(Nodes.if_r_minus_index, Nodes.if_i_minus_index))]) > self.NR_tolerance):
                        print('Stationarity condition (dL/dif) violated')
                        violation = 1
                    if np.any(prim_feas_flag) == True:
                        print('Primal feasibility (inequality) violated at nodes', np.where(prim_feas_flag == True))
                        violation = 1
                    if np.any(dual_feas_flag) == True:
                        print('Dual feasibility (inequality) violated at nodes', np.where(dual_feas_flag == True))
                        violation = 1
                else:
                    if np.any(abs(res_eqn[Nodes.if_index]) > self.NR_tolerance):
                        print('Stationarity condition (dL/dif) violated')

            elif self.source_type == 'PQ' or self.source_type == 'Q':
                if np.any(abs(res_eqn[Nodes.infeas_imag_var_index]) > self.NR_tolerance):
                    print('Stationary constraint (slack Q) violated')
                    violation = 1

                if self.source_type == 'PQ':
                    if np.any(abs(res_eqn[Nodes.infeas_real_var_index]) > self.NR_tolerance):
                        print('Stationary constraint (slack P) violated')
                        violation = 1

                if self.obj == 'L1':
                    if np.any(V_out[Nodes.mu_index] < 0):
                        print('Dual feasibility (inequality) violation; mu is negative')
                        violation = 1

                    if np.any(abs(res_eqn[Nodes.mu_index]) > self.NR_tolerance):
                        print('Complementary slackness condition violated')
                        violation = 1

                    if np.any(V_out[Nodes.infeas_imag_var_index] < 0):
                        print('Primal feasibility (inequality) violated; an infeas Q is negative')
                        violation = 1

                    if self.source_type == 'PQ':
                        if np.any(V_out[Nodes.infeas_real_var_index] < 0):
                            print('Primal feasibility (inequality) violated; an infeas P is negative')
                            violation = 1
            elif self.source_type == 'GB' or self.source_type == 'B':
                print('line 191 InfeasibilityAnalysis - need to build out residuals')

            if self.voltage_unbalance_eqn_enabled:
                if np.any(V_out[Nodes.voltage_unbalance_mu_nodes]) <= 0:
                    print('Voltage Unbalance dual feasibility (inequality) violated')

        if violation == 0:
            print("No violations from residuals")

        if violation == 1:
            list_of_nodes = []
            violated_indices = np.where(abs(res_eqn) >= self.NR_tolerance)[0]

            for i in range(len(self.node)):
                node_set = np.array(self.node[i].node_set)
                if any(node_num in violated_indices for node_num in node_set):
                    list_of_nodes.append(self.node[i].name)
        return res_eqn

    def run_NR_loop(self, V, Ylin, idx_Ynlin, idx_Ylin_H, Jlin, sizeY, errMax_Store, rng, stamped_ground):
        err_max = 1.1 * self.NR_tolerance
        self.innerLoopCount = 0

        # Flags
        flag_isnan = False
        flag_maxiter = False
        flag_tx = False

        while err_max > self.NR_tolerance:
            # Increase iteration count
            self.innerLoopCount += 1
            iteration_count = next(self._iteration_count)

            # if turnoffVlim is True:
            #     voltage_limiting = False

            if iteration_count > self.max_iter:
                flag_maxiter = True
                break

            """ ===========================STAMP Non-Linear ==================="""
            # Stamp Non-Linear Elements
            (Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin,
                idx_Jnlin) = stamp_nonlinear(self.node_key, self.node, self.regulator, self.load, self.triplex_load, self.ibdg,
                                            self.voltage_bounds, self.curr_meas, self.features, V, idx_Ynlin,
                                            self.load_factor, self.homotopy_enabled, self.homotopy, True)

            if self.stamp_dual and (self.obj == 'L1' and self.source_type == 'current'):
                for ele in range(len(self.infeasibility_sources)):
                    infeas_connected_node = self.node_key[self.infeasibility_sources[ele].connected_node_name]
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.infeasibility_sources[ele].stamp_nonlinear_infeasibility(
                                                                                                                self.node[infeas_connected_node], V, self.cs_eps,
                                                                                                                Ynlin_val, Ynlin_row, Ynlin_col,
                                                                                                                idx_Ynlin, Jnlin_val, Jnlin_row,
                                                                                                                idx_Jnlin)

            if self.stamp_dual and (self.source_type == 'PQ' or self.source_type == 'Q'):
                for ele in range(len(self.infeasibility_sources)):
                    infeas_connected_node = self.node_key[self.infeasibility_sources[ele].connected_node_name]
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.infeasibility_sources[ele].stamp_nonlinear_infeasibility_PQ(
                                                                                                                self.node[infeas_connected_node], V,
                                                                                                                Ynlin_val, Ynlin_row, Ynlin_col,
                                                                                                                idx_Ynlin, Jnlin_val, Jnlin_row,
                                                                                                                idx_Jnlin, self.cs_eps)

            if self.stamp_dual and (self.source_type == 'GB' or self.source_type == 'B'):
                for ele in range(len(self.infeasibility_sources)):
                    infeas_connected_node = self.node_key[self.infeasibility_sources[ele].connected_node_name]
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.infeasibility_sources[ele].stamp_nonlinear_infeasibility_GB(self.infeasibility_sources[ele],
                                                                                                                self.node[infeas_connected_node], V,
                                                                                                                Ynlin_val, Ynlin_row, Ynlin_col,
                                                                                                                idx_Ynlin, Jnlin_val, Jnlin_row,
                                                                                                                idx_Jnlin)

            # BATTERY SOURCE STAMPING

            if (self.stamp_dual):
                for ele in range(len(self.battery_sources)):
                    battery_connected_node = self.node_key[self.battery_sources[ele].connected_node_name]
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.battery_sources[ele].stamp_nonlinear_infeasibility(
                                                                                                                self.node[battery_connected_node], V,
                                                                                                                Ynlin_val, Ynlin_row, Ynlin_col,
                                                                                                                idx_Ynlin, Jnlin_val, Jnlin_row,
                                                                                                                idx_Jnlin, self.cs_eps)


            # VOLTAGE UNBALANCE
            if self.voltage_unbalance_eqn_enabled:
                for ele in self.node:
                    if ele.phases == 7 or ele.phases == 15:
                        Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.voltage_unbalance.stamp_voltage_unbalance_constraints(ele, V, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin)

            # Create Y and J nonlin matrix
            Jnlin_col = np.zeros((len(Jnlin_row)), dtype=np.int64)

            Ynlin, Jnlin = matrix_conversion(Ynlin_val, Ynlin_row, Ynlin_col,
                                                Jnlin_row, Jnlin_col, Jnlin_val,
                                                sizeY)
            """===============================End==============================="""
            """---------------------Homotopy Stepping---------------------------"""
            # Create linear homotopy vectors

            Ylin_val_H = np.zeros(2*idx_Ylin_H, dtype = float)
            Ylin_col_H = np.zeros(2*idx_Ylin_H, dtype = int)
            Ylin_row_H = np.zeros(2*idx_Ylin_H, dtype = int)

            Jlin_val_H = np.zeros(2*idx_Ylin_H, dtype = float)
            Jlin_row_H = np.zeros(2*idx_Ylin_H, dtype = int)
            Jlin_col_H = np.zeros(2*idx_Ylin_H, dtype = int)

            # Stamp homotopy
            if self.homotopy_enabled and self.homotopy.h_factor != 0:
                Ylin_row_H, Ylin_col_H, Ylin_val_H, idx_Ylin_H, stamped_ground = self.homotopy.stamp_tx_homotopy(
                                            self.node, self.node_key, self.line_oh, self.line_ug,
                                            self.line_tplx, self.xfmr, stamped_ground, Ylin_val_H,
                                            Ylin_row_H, Ylin_col_H, self.obj, Jlin_val_H, Jlin_row_H, V)

            #Ylin_H = csr_matrix((Ylin_val_H, (Ylin_row_H, Ylin_col_H)), shape=(sizeY, sizeY), dtype=np.float64)
            Ylin_H, Jlin_H = matrix_conversion(Ylin_val_H, Ylin_row_H, Ylin_col_H,
                                                Jlin_row_H, Jlin_col_H, Jlin_val_H,
                                                sizeY)

            """====================Combine The linear and Non-linear parts======"""
            Y = Ylin + Ynlin + Ylin_H
            J = Jlin + Jnlin + Jlin_H

            # Reformat the Y matrix to not include rows and cols of zeros
            (Y_red, J_red, sol_index) = check_zero_rows(Y, J)

            """ ==================Solve the system and calculate the error======"""
            Vsol_red = spsolve(Y_red, J_red)
            Vsol_red = np.reshape(Vsol_red, (Vsol_red.shape[0], 1))
            # Map the Vsol solution to the original V vector
            Vsol = np.zeros((sizeY, 1), dtype=float)
            Vsol[sol_index] = Vsol_red

            # Calculate the error

            Verr = V#[Nodes.voltage_index]

            Vsol_err = Vsol#[Nodes.voltage_index]
            err = Verr - Vsol_err
            err_max = np.amax(abs(err))
            errMax_Store.append(err_max)

            # Figure out what node has the error
            argmax_ind = np.argmax(abs(err))
            err_max_node = [
                ele.name for ele in self.node if argmax_ind in ele.node_set
            ]
            if err_max_node:
                print(
                    colored(
                        'Maximum error from iteration %s is %f at node %s' %
                        (iteration_count, err_max, err_max_node[0]), 'red'))
            else:
                print(
                    colored('Maximum error from iteration %s is %f' % (iteration_count, err_max),
                            'red'))
            if math.isnan(err_max):
                flag_isnan = True
                break

            # if np.any(Vsol[Nodes.voltage_unbalance_tracking_nodes] > (self.voltage_unbalance.upper_bound**2/100**2)):
            #     print('Vunb violates constraint')
            """--------------------  Tx Steppping ------------------------------"""
            # Run tx_stepping
            if self.homotopy_enabled:
                Vsol = self.homotopy.run_tx_stepping(err_max, Vsol)

                # If homotopy failed, break the loop
                if self.homotopy.tx_stepping_failed == True:
                    flag_tx = True
                    break

                if self.homotopy.perturb_conditions == True:
                    random_ones = rng.random(len(Vsol))
                    random_ones[random_ones > 0.5] = 1
                    random_ones[random_ones <= 0.5] = -1
                    random_ones = random_ones.reshape((len(random_ones),1))

                    if self.obj == 'L1':
                        # temp_indices = np.where(Vsol[self.node[0].if_index] > 1e-6)
                        # new_Vinit = Vsol + random_ones*0.15*Vsol
                        #new_Vinit[self.node[0].if_index] += 1
                        random_infeas_ones = rng.random(len(Vsol[self.node[0].if_index]))
                        random_infeas_ones[random_infeas_ones > 0.25] = 10
                        random_infeas_ones[random_infeas_ones <= 0.25] = 0
                        new_Vinit = Vsol # + random_ones * 0.15 * Vsol
                        new_Vinit[self.node[0].if_index] += random_infeas_ones.reshape((len(random_infeas_ones),1)) #np.reshape(random_infeas_ones, (len(random_infeas_ones), 1))
                    elif self.obj == 2:
                        new_Vinit = Vsol + random_ones*0.15*Vsol + 1e-3
                    else:
                        new_Vinit = Vsol + random_ones*0.15*Vsol
                    Vsol = new_Vinit

                    self.homotopy.perturb_conditions = False

            """--------------------  Voltage Limiting -------------------------"""
            if self.voltage_limiting_enabled and err_max > self.NR_tolerance:
                V = self.voltage_limiting.apply_voltage_limiting(Vsol, V)

            """--------------------  Voltage Limiting -------------------------"""
            if self.diode_limiting_enabled and err_max > self.NR_tolerance:
                if (self.obj == 'L1' and self.source_type == 'current'):
                    Vsol = self.diode_limiting_infeasibility.apply_diode_limiting(Vsol, V)
                    if np.size(np.where(V[Nodes.mu_index] < 0)[0]) != 0:
                        print('Mus are negative')
                elif (self.obj == 'L1' and self.source_type == 'Q'):
                    Vsol = self.diode_limiting_infeasibility_q.apply_diode_limiting(Vsol, V)
                    if np.size(np.where(V[Nodes.mu_index] < 0)[0]) != 0:
                        print('Mus are negative')
                elif (self.obj == 'L1' and self.source_type == 'PQ'):
                    Vsol = self.diode_limiting_infeasibility_q.apply_diode_limiting(Vsol, V)
                    Vsol = self.diode_limiting_infeasibility_p.apply_diode_limiting(Vsol,V)
                    if np.size(np.where(V[Nodes.mu_index] < 0)[0]) != 0:
                        print('Mus are negative')

                if self.voltage_unbalance_eqn_enabled:
                    Vsol = self.diode_limiting_unbalance.apply_diode_limiting_upper_bound_only(Vsol, V)
                    if np.size(np.where(V[Nodes.voltage_unbalance_mu_nodes] < 0)[0]) != 0:
                        print('Mus are negative')
                    if np.any(Vsol[Nodes.voltage_unbalance_tracking_nodes] > (self.voltage_unbalance.upper_bound**2/100**2)):
                        print('Vunb violates constraint')

            if self.voltage_bounds and err_max > self.NR_tolerance:
                V_limited = operational_limits_diode_limiting(Vsol, V, self.voltage_bound_settings["cs eps tol"], 'bus_voltages')
                Vsol = np.copy(V_limited)


            """--------------------  Battery Power & SOC Limiting -------------------------"""

            for bat in self.battery_sources:
                # Perform Upper and Lower Limiting for Power Distribution
                mu_p_upper = bat.mu_index_upper
                P_index = bat.P_plus_nodes + bat.P_minus_nodes
                Vsol = operational_limits_diode_limiting(Vsol, V, cs_eps = 1e-5, type = 'battery', d=0.95, normalize=False,
                                                         bat = bat)

                # Perform Upper and Lower Limiting for SOC Distribution
                mu_Bt_upper = bat.mu_index_Bt_upper
                Bt_index = bat.Bt_nodes




            V = np.copy(Vsol)
            """--------------------  Residuals -------------------------"""
            # if iteration_count % 10 == 0:
            #	res_eqn = self.calc_residuals(V, Ylin, Jlin, self.node, self.node_key, self.load, self.triplex_load, self.load_factor, True, self.obj, self.cs_eps)

            """=========================End====================================="""
            if not self.homotopy_enabled and err_max < self.NR_tolerance:
                break

        return Vsol, err_max, flag_isnan, flag_tx, flag_maxiter

    def run_infeasibility_analysis(self):
        returncode = {'success': 0, 'fail': 1, 'error': 2}
        flag_tx = False

        simulation_start_time = time.perf_counter()
        # Size of the solution vector
        sizeY = int(repr(self.node_index_)[6:-1])
        print("size Y is ", sizeY)
        # Initialize V vector
        Vinit = np.zeros([sizeY, 1], dtype=np.float64)
        """------------------------------End----------------------------------------"""
        """-----------------------Initialize V and Q Vector-------------------------"""

        Vinit = initialize(self.node, self.node_key, self.regulator, self.ibdg, self.battery_sources, self.voltage_bounds, Vinit, self.load_factor, self.stamp_dual, self.obj, self.source_type, self.stamp_slack_flag, self.stamp_tplx_infeas_sources_flag, self.stamp_neutral_infeas_source_flag)
        
        #
        if self.use_Vinit:
            Vinit = self.Vinit
        if self.warm_start:
            if self.settings["Warm Start"]["from external simulation"]:
                Vinit = read_open_dss_results(self.settings["Warm Start"]["solution file path"], Vinit, self.node, self.node_key)
                # To-do read in current file and intialize all of the currents in the system
                Vinit[self.node[0].slack_nodeA_Vr] = 267.416*math.cos(math.degrees(-172.1))
                Vinit[self.node[0].slack_nodeA_Vi] = 267.416*math.sin(math.degrees(-172.1))
                Vinit[self.node[0].slack_nodeB_Vr] = 281.359*math.cos(math.degrees(65.42))
                Vinit[self.node[0].slack_nodeB_Vi] = 281.359*math.sin(math.degrees(65.42))
                Vinit[self.node[0].slack_nodeC_Vr] = 267.742*math.cos(math.degrees(-52.52))
                Vinit[self.node[0].slack_nodeC_Vi] = 267.742*math.sin(math.degrees(-52.52))
            else:
                warm_start_settings = self.settings["Warm Start"]
                Vinit = initialize_warm_start(Vinit, self.case_name, self.load_factor, self.node, warm_start_settings)

        if self.voltage_unbalance_eqn_enabled:
            for ele in self.node:
                Vinit = self.voltage_unbalance.initialize_voltage_unbalance_variables(ele, Vinit)
        # Initializing initial size
        if sizeY < 1000:
            idx_Ylin = 100 * sizeY
            idx_Ynlin = 100 * sizeY
            idx_Ylin_H = 100 * sizeY
        else:
            idx_Ylin = 10 * sizeY
            idx_Ynlin = 10 * sizeY
            idx_Ylin_H = 10 * sizeY

        print(
            '======================= Variable Initialization Complete =============================='
        )
        """------------------------------End----------------------------------------"""

        """ ===========================STAMP Linear ================================"""
        stamped_ground = set()

        # Stamp Linear Elements
        (Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Ylin,
         idx_Jlin, stamped_ground) = stamp_linear(self.node_key, self.slack, self.line_oh,
                                self.line_ug, self.line_tplx,
                                self.node, self.xfmr, self.regulator, self.switch, self.fuse, self.capacitor,
                                self.reactor, self.curr_meas, self.features, stamped_ground, idx_Ylin, Vinit,
                                self.stamp_dual)

        if self.stamp_dual and self.source_type == 'current':
            for ele in range(len(self.infeasibility_sources)):
                infeas_connected_node = self.node_key[self.infeasibility_sources[ele].connected_node_name]
                Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Ylin, idx_Jlin = self.infeasibility_sources[ele].stamp_linear_infeasibility(
                                                                                                            self.node[infeas_connected_node], Ylin_val,
                                                                                                            Ylin_row, Ylin_col,
                                                                                                            idx_Ylin, Jlin_val,
                                                                                                            Jlin_row, idx_Jlin)

        Jlin_col = np.zeros((len(Jlin_row)), dtype=np.int64)
        Ylin, Jlin = matrix_conversion(Ylin_val, Ylin_row, Ylin_col, Jlin_row,
                                       Jlin_col, Jlin_val, sizeY)
        res_eqn = np.zeros_like(Vinit)

        # Raw Variables
        errMax_Store = [1e-1]
        num_outer = 0
        self._iteration_count = count(0)

        # Flags
        flag_maxiter = False
        flag_isnan = False

        # Map System
       #if self.draw_network_maps:
        #    Map.make_KVA_map(self.node, self.line_oh, self.line_ug, self.line_tplx, self.xfmr, self.regulator, self.fuse, self.switch)
        #    Map.make_map(self.node, self.line_oh, self.line_ug, self.line_tplx, self.xfmr, self.regulator, self.fuse, self.switch)

        if self.diode_limiting_enabled:
            if (self.obj == 'L1' and self.source_type == 'current'):
                if_index = Nodes.if_index
                mu_index = Nodes.mu_index
                diode_limiting_dict_infeas = {'xmin': np.zeros((np.size(if_index,0),1)), 'x_index': if_index, 'mu_min_index': mu_index, 'cs_tol': self.cs_eps}
                self.diode_limiting_infeasibility = Limiting('diode', diode_limiting_dict_infeas)

            if (self.obj == 'L1' and self.source_type == 'Q'):
                Q_index = Nodes.infeas_imag_var_index
                mu_q_index = Nodes.mu_i_index
                diode_limiting_dict_infeas = {'xmin': np.zeros((np.size(Q_index,0),1)), 'x_index': Q_index, 'mu_min_index': mu_q_index, 'cs_tol': self.cs_eps}
                self.diode_limiting_infeasibility_q = Limiting('diode', diode_limiting_dict_infeas)

            if (self.obj == 'L1' and self.source_type == 'PQ'):
                Q_index = Nodes.infeas_imag_var_index
                mu_q_index = Nodes.mu_i_index
                diode_limiting_dict_infeas_q = {'xmin': np.zeros((np.size(Q_index,0),1)), 'x_index': Q_index, 'mu_min_index': mu_q_index, 'cs_tol': self.cs_eps}
                self.diode_limiting_infeasibility_q = Limiting('diode', diode_limiting_dict_infeas_q)

                # mu is constrained to be larger than zero
                P_index = Nodes.infeas_real_var_index
                mu_p_index = Nodes.mu_r_index
                diode_limiting_dict_infeas_p = {'xmin': np.zeros((np.size(P_index,0),1)), 'x_index': P_index, 'mu_min_index': mu_p_index, 'cs_tol': self.cs_eps}
                self.diode_limiting_infeasibility_p = Limiting('diode', diode_limiting_dict_infeas_p)

            if self.voltage_unbalance_eqn_enabled:
                v_unb_index = Nodes.voltage_unbalance_tracking_nodes
                mu_index = Nodes.voltage_unbalance_mu_nodes
                max_unbalance = np.ones((np.size(v_unb_index,0),1))*(self.voltage_unbalance.upper_bound**2/100**2)
                diode_limiting_dict_unbalance = {'xmax': max_unbalance, 'x_index': v_unb_index, 'mu_max_index': mu_index, 'cs_tol': self.cs_eps}
                self.diode_limiting_unbalance = Limiting('diode', diode_limiting_dict_unbalance)

            # lower diode limiting for batteries
            """
            for bat in self.battery_sources:
                P_index = bat.P_plus_nodes + bat.P_minus_nodes
                mu_p_index = bat.mu_index
                diode_limiting_battery_p = {'xmin': np.zeros((np.size(P_index,0),1)), 'x_index': P_index, 'mu_min_index': mu_p_index, 'cs_tol': self.cs_eps}
			"""


        if self.homotopy_enabled:
            self.homotopy.Vinit = Vinit
        V = np.copy(Vinit)
        rng = np.random.default_rng(0)

        while True:
            """ Variable for the While Loop"""
            num_outer += 1
            err_max = 1.1 * self.NR_tolerance
            res_eqn = self.calc_residuals(Vinit, Ylin, Jlin, self.node_key)


            Vsol, err_max, flag_isnan, flag_tx, flag_maxiter = self.run_NR_loop(V, Ylin, idx_Ynlin, idx_Ylin_H, Jlin, sizeY, errMax_Store, rng, stamped_ground)
            
			# display battery outputs for testing
            if (self.battery_sources):
                for i,bat in enumerate(self.battery_sources):
                    if (bat.verbose):
                        bat.show_power_output(Vsol)


            """##############OUTER LOOP OF POWERFLOW OPTIONS#################"""
            """Mininum power_step or minimum iter fatal stop"""

            if flag_maxiter:
                logging.debug(
                    "Simulation has failed. Maximum iteration count reached.")
                exit()
                break

            if flag_tx:
                logging.debug(
                    "Simulation has failed with homotopy.")
                exit()
                break

            if flag_isnan:
                logging.debug("Simulation has failed.")
                exit()
                break

            if self.homotopy_enabled:
                if err_max < self.NR_tolerance and self.homotopy.tx_stepping_successful == True:
                    break
            else:
                if err_max < self.NR_tolerance:
                    break

        V = Vsol

        simulation_end_time = time.perf_counter()
        sim_total_time = simulation_end_time - simulation_start_time
        self.simulation_stats.append(sim_total_time)
        self.Vsol = V
        logging.info(colored('Simulation time is %f seconds' % (sim_total_time), 'blue'))

        if flag_maxiter:
            logging.error('Case did not converge - Maximum Iteration Reached')
            return returncode["fail"]
        elif flag_tx:
            logging.error('Case did not converge - Homotopy failed')
            return returncode["fail"]
        else:
            iteration_count = next(self._iteration_count)
            logging.info(
                colored(
                    'Case converged in %d iterations - SUCCESS!' % iteration_count,
                    'green'))

            if self.save_output_conditions_flag:
                if self.obj == 'L2':
                    self.save_output_conditions("_L2node_", "_L2_V_", Vsol)

                if self.obj == 'L1':
                    self.save_output_conditions("_L1node_", "_L1_V_", Vsol)

                if self.obj == None:
                    self.save_output_conditions("_pf_node_", "_pf_V_", Vsol)

            self.res_eqn = self.calc_residuals(Vsol, Ylin, Jlin, self.node_key)

            if self.normalize_if:
                normalized_if, max_if = normalize_if(self.node, Vsol, self.obj)

            if self.report_if and self.stamp_dual:
                dual_info_ABC, dual_info_tplx = dual_analysis(Vsol, self.node, self.obj, self.stamp_slack_flag, self.stamp_tplx_infeas_sources_flag, self.source_type)
            else:
                dual_info_ABC = None
                dual_info_tplx = None

            calc_voltage_violations(Vsol, self.node, self.vmax_pu, self.vmin_pu)

            if self.draw_infeasibility_maps:
                Map.map_infeasibility(self.node, self.xfmr, self.regulator, self.switch, self.load,
                                        self.line_oh, self.line_ug, self.line_tplx, self, self.fuse, Vsol)
                if self.normalize_if:
                    Map.map_normal_infeasibility(self.node, self.xfmr, self.regulator, self.switch, self.load,
                                            self.line_oh, self.line_ug, self.line_tplx, self, self.fuse, Vsol, max_if)
            if self.save_results:
                self.write_results(V, iteration_count, self.simulation_stats, dual_info_ABC, dual_info_tplx)

            return returncode["success"]
