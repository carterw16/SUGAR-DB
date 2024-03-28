"""Three-phase power flow analysis

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-10-2017
  Updated Date: 10-14-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Validated

"""

import logging
import numpy as np
import time
from lib.stamp_nonlinear import stamp_nonlinear
from classes.CurrentMeasurements import create_current_meas
from classes.Homotopy import Homotopy
from lib.matrix_conversion import matrix_conversion
from lib.initialize import initialize
from lib.stamp_linear import stamp_linear
from termcolor import colored
from itertools import count
from scipy.sparse import csr_matrix
from lib.check_zero_rows import check_zero_rows
from lib.check_error import check_error
from scipy.sparse.linalg import spsolve
from classes.Nodes import Nodes
from lib.apply_voltage_limiting import apply_voltage_limiting
from lib.save_output import save_output
from lib.assign_nodes import assign_nodes


class PowerFlow:
    def __init__(self, case_name, casedata, features, settings, path_to_output, node_key):

        self.path_to_output = path_to_output
        self.case_name = case_name
        self.simulation_stats = casedata.stats
        self.save_results = settings["Save File"]

        """----------------------Assign nodes to objects-------------------------"""
        (node_index_, casedata.node,
         casedata.regulator, casedata.xfmr, casedata.switch,
         casedata.fuse, casedata.reactor, casedata.ibdg) = assign_nodes(node_key, casedata.node, casedata.regulator,
                                                                        casedata.xfmr, casedata.switch,
                                                                        casedata.fuse, casedata.reactor,
                                                                        casedata.ibdg)

        self.nodes = casedata.node
        self.loads = casedata.load
        self.slacks = casedata.slack
        self.triplex_loads = casedata.triplex_load
        self.transformers = casedata.xfmr
        self.oh_lines = casedata.ohline
        self.ug_lines = casedata.ugline
        self.triplex_lines = casedata.triplex_line
        self.capacitors = casedata.capacitor
        self.regulators = casedata.regulator
        self.switches = casedata.switch
        self.fuses = casedata.fuse
        self.reactors = casedata.reactor
        self.ibdgs = casedata.ibdg

        self.node_index_ = node_index_
        self.node_key = node_key

        self.features = features
        self.settings = settings

        self.lf = settings["Load Factor"]
        self.tol = settings["Tolerance"]
        self.max_iter = settings["Max Iters"]
        self.num_outer = 0
        self.voltage_limiting = settings["Voltage Limiting"]
        self.flag_maxiter = False
        self.inner_loop_complete = True

        # # Homotopy # #
        self.flag_tx = False
        self.homotopy_enabled = settings["Homotopy"]
        g_homotopy = settings["G_homotopy"]
        b_homotopy = settings["B_homotopy"]
        h_factor_init = 1
        homotopy_tolerance = 1e-4
        h_step_size_init = 0.1

        self.homotopy = Homotopy(homotopy_tolerance, g_homotopy, b_homotopy,
                                 h_factor_init, h_step_size_init)

        # # Current Measurements # #
        # Create ammeters and places them in a current measurement list curr_meas if enabled
        # The ammeters measure overhead, underground, and triplex line currents
        enable_CM = self.features['Current Meas']
        self.curr_meas = []
        if enable_CM:
            self.activate_current_measures(enable_CM, self.node_index_, self.node_key)

        # # Create Solution Vector # #
        # Size of the solution vector
        self.size_Y = int(repr(self.node_index_)[6:-1])
        logging.debug("size Y is %d" % self.size_Y)
        # Initialize V vector
        self.V_init = np.zeros([self.size_Y, 1], dtype=np.float)
        self.V_init = initialize(self.node_key, self.nodes, self.regulators, self.ibdgs, self.V_init, self.lf)

    def activate_current_measures(self, enable_CM, node_index_, node_key):
        self.curr_meas = create_current_meas(enable_CM, self.oh_lines, self.ug_lines,
                                             self.triplex_lines)

        for ele in self.curr_meas:
            node_index_ = ele.assign_nodes(node_key, node_index_)

        self.node_index_ = node_index_
    
    def calc_residuals(V_out, Ylin, Jlin, nodes, node_key, load, triplex_load, lf, stamp_dual, obj, epsilon):
        '''

        Parameters
        ----------
        V_out : converged solution from Newton Raphson iterations
        Ylin : linear Y matrix
        Jlin : linear J vector
        nodes : list of node elements
        load : list of load elements
        triplex_load : list of triplex load elements
        lf : scalar load factor
        stamp_dual : boolean variable for whether or not to incorporate dual equations
        obj : scalar of value 1 or 2; built out for either L1 norm or L2 norm infeasibility analysis
        epsilon : scalar value (should be small) used as relaxation in complementary slackness

        Returns
        -------
        res_eqn : vector the size of V_out; the result of every original KKT condition 
            evaluated at the optimal values from V; individual elements should equal 0

        '''
    # Initializing the vector which houses the residual equations
        res_eqn = np.zeros_like(V_out)
        
        violation = 0
        # Primal Feasibility (Equality Constraints) [some elements in this equation are nonlinear, so these components
        # are directly added to res_eqn for load and triplex load]
        # eq constr = 0 (located in V nodes) 
        for ele in range(len(load)):
            res_eqn = load[ele].calc_residual(V_out, node_key, nodes, res_eqn, lf)


        for ele in range(len(triplex_load)):
            res_eqn = triplex_load[ele].calc_residual(V_out, node_key, nodes, res_eqn, lf)
        
        prim_feas_flag = np.empty(len(nodes))

        # If res_eqn = 0, then stationarity, primal feasibility (equality), and complementary slackness are met
        res_eqn += (Ylin@V_out - Jlin)

        if any(abs(res_eqn[np.append(Nodes.Vr_index, Nodes.Vi_index)]) > 1e-6):
            print('Primal feasibility (equality) condition violated')
            violation = 1
            
        if violation == 0:
            print("No violations from residuals")
            
        if violation == 1:
            list_of_nodes = []
            violated_indices = np.where(abs(res_eqn) >= 1e-8)[0]
            
            for i in range(len(nodes)):
                node_set = np.array(nodes[i].node_set)
                if any(node_num in violated_indices for node_num in node_set):
                    list_of_nodes.append(nodes[i].name)
        return res_eqn

    def run_pf(self):
        simulation_start_time = time.perf_counter()

        # Declare return codes
        return_code = {'success': 0, 'fail': 1, 'error': 2}

        # Voltage Limiting
        max_step__ = 0.1
        min_step__ = -0.1
        Vr_limit__ = 2
        Vi_limit__ = 2

        # Raw Variables
        Vr_norm = None
        Vi_norm = None
        err_max_store = [1e-1]
        lfStore = []

        # Initializing initial size
        if self.size_Y < 1000:
            idx_Ylin = 100 * self.size_Y
            idx_Ynlin = 100 * self.size_Y
            idx_Ylin_H_init = 100 * self.size_Y
        else:
            idx_Ylin = 50 * self.size_Y
            idx_Ynlin = 50 * self.size_Y
            idx_Ylin_H_init = 50 * self.size_Y

        logging.debug(
            '======================= Variable Initialization Complete =============================='
        )
        """------------------------------End----------------------------------------"""

        """ ===========================STAMP Linear ================================"""
        # Track Stamped Ground Nodes
        stamped_ground = set()
        # Stamp Linear Elements
        (Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Ylin, idx_Jlin,
         stamped_ground) = stamp_linear(self.node_key, self.slacks, self.oh_lines, self.ug_lines,
                                        self.triplex_lines, self.nodes, self.V_init, self.transformers,
                                        self.regulators,
                                        self.switches, self.fuses, self.capacitors, self.reactors, self.curr_meas,
                                        self.features, stamped_ground,
                                        idx_Ylin)

        Jlin_col = np.zeros((len(Jlin_row)), dtype=np.int64)
        Ylin, Jlin = matrix_conversion(Ylin_val, Ylin_row, Ylin_col, Jlin_row,
                                       Jlin_col, Jlin_val, self.size_Y)

        V = np.copy(self.V_init)

        _iteration_count = count(0)

        logging.info(
            colored(
                '======================= Start SUGAR3 PF Solver ==============================',
                'white'))

        while True:

            """ Variable for the While Loop"""
            self.num_outer += 1
            err_max = 1.1 * self.tol
            inner_loop_count = 0

            while err_max > self.tol:

                # Increase iteration count
                inner_loop_count += 1
                iteration_count = next(_iteration_count)

                if iteration_count > self.max_iter:
                    flag_maxiter = True
                    break

                """ ===========================STAMP Non-Linear ==================="""
                # Stamp Non-Linear Elements
                (Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin,
                 idx_Jnlin) = stamp_nonlinear(self.node_key, self.nodes, self.regulators, self.loads,
                                              self.triplex_loads, self.ibdgs, self.curr_meas, self.features,
                                              V, idx_Ynlin, self.lf, self.homotopy_enabled, self.homotopy)

                # Create Y and J nonlin matrix
                Jnlin_col = np.zeros((len(Jnlin_row)), dtype=np.int64)

                Ynlin, Jnlin = matrix_conversion(Ynlin_val, Ynlin_row, Ynlin_col,
                                                 Jnlin_row, Jnlin_col, Jnlin_val,
                                                 self.size_Y)
                """===============================End==============================="""
                """---------------------Homotopy Stepping---------------------------"""
                # Create linear homotopy vectors
                Ylin_val_H = np.zeros(2 * idx_Ylin_H_init, dtype=float)
                Ylin_col_H = np.zeros(2 * idx_Ylin_H_init, dtype=int)
                Ylin_row_H = np.zeros(2 * idx_Ylin_H_init, dtype=int)

                # Stamp homotopy
                if self.homotopy_enabled and self.homotopy.h_factor != 0:
                    Ylin_row_H, Ylin_col_H, Ylin_val_H, idx_Ylin_H, stamped_ground = self.homotopy.stamp_tx_homotopy(
                        self.node_key, self.nodes, self.oh_lines, self.ug_lines, self.triplex_lines, self.transformers,
                        stamped_ground, Ylin_val_H, Ylin_row_H, Ylin_col_H)

                Ylin_H = csr_matrix((Ylin_val_H, (Ylin_row_H, Ylin_col_H)),
                                    shape=(self.size_Y, self.size_Y),
                                    dtype=np.float64)

                """====================Combine The linear and Non-linear parts======"""
                Y = Ylin + Ynlin + Ylin_H
                J = Jlin + Jnlin

                # Reformat the Y matrix to not include rows and cols of zeros
                (Y_red, J_red, sol_index) = check_zero_rows(Y, J)
                """ ==================Solve the system and calculate the error======"""
                Vsol_red = spsolve(Y_red, J_red, use_umfpack=True)

                Vsol_red = np.reshape(Vsol_red, (Vsol_red.shape[0], 1))
                # Map the V_sol solution to the original V vector
                V_sol = np.zeros((self.size_Y, 1), dtype=float)
                V_sol[sol_index] = Vsol_red

                # # Calculate the iteration error and pinpoint the location # #
                err_nodes = []
                voltage_index = Nodes.voltage_index
                tap_index = Nodes.Tap_index
                Q_index = Nodes.Q_index

                err_nodes.extend(voltage_index)
                err_nodes.extend(tap_index)
                err_nodes.extend(Q_index)

                err_max, err_max_store, err_max_node = check_error(err_nodes, V, V_sol, err_max_store, self.nodes,
                                                                   self.regulators,
                                                                   self.ibdgs)
                logging.info(
                    colored(
                        'Maximum error from this iteration is %f at node %s of type %s' %
                        (err_max, err_max_node[0][0], err_max_node[0][1]), 'white'))

                if np.isnan(err_max):
                    break

                lfStore.append(self.lf)

                """--------------------  Tx Steppping ------------------------------"""
                # Run tx_stepping
                if self.homotopy_enabled:
                    self.homotopy.run_tx_stepping(err_max, V_sol)

                    # If homotopy failed, break the loop
                    if self.homotopy.tx_stepping_failed == True:
                        flag_tx = True
                        break

                """--------------------  Voltage Limiting -------------------------"""
                if self.voltage_limiting and err_max > self.tol:
                    apply_voltage_limiting(V_sol, V, Nodes, max_step__, min_step__,
                                           Vr_norm, Vi_norm, Vr_limit__, Vi_limit__)
                else:
                    V = np.copy(V_sol)
                """=========================End====================================="""

                """=========================End====================================="""

                if err_max < self.tol:
                    break

            """##############OUTER LOOP OF POWERFLOW OPTIONS#################"""
            """Mininum power_step or minimum iter fatal stop"""
            if self.flag_maxiter:
                logging.info(
                    "Simulation has failed. Maximum iteration count reached.")
                exit()
                break
            if self.flag_tx:
                logging.error('Case did not converge - Homotopy did not converge')
                exit()
                break

            if self.inner_loop_complete:
                self.inner_loop_complete = False

            # Finish The Simulation if:
            # 1. Max error less than tolerance
            if self.homotopy_enabled:
                if err_max < self.tol and self.homotopy.tx_stepping_successful:
                    break
            else:
                if err_max < self.tol:
                    break

        # Print End time of the simulation
        simulation_end_time = time.perf_counter()
        sim_total_time = simulation_end_time - simulation_start_time
        self.simulation_stats.append(sim_total_time)
        logging.info(colored('Simulation time is %f seconds' % (sim_total_time), 'blue'))

        if self.flag_maxiter:
            logging.error('Case did not converge - Maximum Iteration Reached')
            return return_code['fail']
        else:
            iteration_count = next(_iteration_count)
            logging.info(
                colored(
                    'Case converged in %d iterations - SUCCESS!' % iteration_count,
                    'green'))
            if self.save_results:
                self.write_results(V, iteration_count, self.simulation_stats)

        return return_code['success']

    def write_results(self, V, iteration_count, simulation_stats):

        for ele in self.nodes:
            ele.print_slack = True
            ele.calcMagAng(V, False)

        if self.regulators:
            for ele in self.regulators:
                ele.update_values(V)

        enable_IBDGs = True if self.ibdgs else False
        if self.ibdgs:
            for ele in self.ibdgs:
                ele.calc_ibdg_outputs(self.nodes, self.node_key, V)

        for ele in self.curr_meas:
            ele.calc_currents(V)

        for ele in self.transformers:
            ele.calc_currents(V)

        for ele in self.fuses:
            ele.calc_currents(V)

        for ele in self.switches:
            ele.calc_currents(V)

        outputs = [
            self.loads, self.ibdgs, self.regulators, self.curr_meas, self.transformers, self.fuses, self.switches,
            self.triplex_loads
        ]

        # Write csv file with results
        save_output(self.path_to_output, self.case_name, self.nodes, outputs, iteration_count,
                    simulation_stats, self.settings, enable_IBDGs, self.lf)
