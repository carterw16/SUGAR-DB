"""
Three-phase battery system classes

Author(s): Meshach Hopkins
Created Date: 02-13-2024
Updated Date: 02-13-2024
Email: mhopkins@andrew.cmu.edu
Status: Development

Parses the power grid elements, initializes their respective classes, and then returns a list for each element.

"""

from types import MethodType
from termcolor import colored
from classes.InfeasibilitySources.base import InfeasibilitySources
from .BatteryStamps import (stampY, stampJ, stamp_nonlinear_infeasibility,
                                    stamp_complementary_slackness, stamp_equality_constraints,
                                    stamp_stationarity_constraints, calculate_partial_derivatives)

import ipdb


class BatterySources(InfeasibilitySources):
    def __init__(self, node, node_index, obj_scalar, source_type, stamp_neutral_infeas_source_flag,
              ID="bat", P_max = 3000000, P_min = 0, Mch = 0.05, Md = 0.05, Bt_prev = 0, Lb_next = 0,
              C_ch = 1, C_d = 1e-10, single_phase = "", verbose = True, bt_mu_tol = 1e-5):
        
        InfeasibilitySources.__init__(self, node)
        self.source_type = source_type
        # the source type for batteries are "P" and "PQ"

        # initialize battery parameters
        self.ID = ID
        self.Lb_next = Lb_next
        self.P_max = P_max
        self.P_min = P_min
        self.Mch = Mch
        self.Md = Md
        # cost function parameters
        self.C_ch = C_ch
        self.C_d = C_d
        # norm type for pch and pd
        self.pch_l2 = True
        self.pd_l2 = False     
        self.single_phase = single_phase
        self.Bt_prev = Bt_prev
        if (self.single_phase == ""):
            if (type(Bt_prev).__name__ == "int" or type(Bt_prev).__name__ == "float" ):
                self.Bt_prev = [Bt_prev, Bt_prev, Bt_prev]
            if (type(Lb_next).__name__ == "int" or type(Lb_next).__name__ == "float" ):
                self.Lb_next = [Lb_next, Lb_next, Lb_next]

        self.verbose = verbose
        self.SOC_max = 100000000 # maximum MWH capacity of battery should be about 100MWH
        self.SOC_min = 0 # default min could change for additional problems

        # adjustable battery mu tolerance variable
        self.bt_mu_tol = 1e-5

        # stamping functions
        self.stamp_neutral = stamp_neutral_infeas_source_flag
        self.stampY = MethodType(stampY, self)
        self.stampJ = MethodType(stampJ, self)
        self.stamp_nonlinear_infeasibility = MethodType(stamp_nonlinear_infeasibility, self)
        self.stamp_complementary_slackness = MethodType(stamp_complementary_slackness, self)
        self.stamp_equality_constraints = MethodType(stamp_equality_constraints, self)
        self.stamp_stationarity_constraints = MethodType(stamp_stationarity_constraints, self)
        self.calculate_partial_derivatives = MethodType(calculate_partial_derivatives, self)

        # assign and update nodes
        node_index = self.assign_nodes(node, node_index)
        self.updated_node_index = node_index

        self.obj_scaling = obj_scalar #1/(node.Vnom**2)


    def initialize(self, Vinit):
        if (self.single_phase == ""):
            Vinit[self.nodeA_p_plus] = 0.0001
            Vinit[self.nodeB_p_plus] = 0.0001
            Vinit[self.nodeC_p_plus] = 0.0001

            Vinit[self.nodeA_p_minus] = 500
            Vinit[self.nodeB_p_minus] = 500
            Vinit[self.nodeC_p_minus] = 500

            Vinit[self.nodeA_dual_ineq_var_p_plus] = 0.0001
            Vinit[self.nodeB_dual_ineq_var_p_plus] = 0.0001
            Vinit[self.nodeC_dual_ineq_var_p_plus] = 0.0001

            Vinit[self.nodeA_dual_ineq_var_p_plus_upper] = 0.0001
            Vinit[self.nodeB_dual_ineq_var_p_plus_upper] = 0.0001
            Vinit[self.nodeC_dual_ineq_var_p_plus_upper] = 0.0001

            Vinit[self.nodeA_dual_ineq_var_p_minus] = 0.0001
            Vinit[self.nodeB_dual_ineq_var_p_minus] = 0.0001
            Vinit[self.nodeC_dual_ineq_var_p_minus] = 0.0001

            Vinit[self.nodeA_dual_ineq_var_p_minus_upper] = 0.0001
            Vinit[self.nodeB_dual_ineq_var_p_minus_upper] = 0.0001
            Vinit[self.nodeC_dual_ineq_var_p_minus_upper] = 0.0001

            if (type(self.Bt_prev).__name__ == 'int' or type(self.Bt_prev).__name__ == 'float'):
                Vinit[self.BtA] = self.Bt_prev
                Vinit[self.BtB] = self.Bt_prev
                Vinit[self.BtC] = self.Bt_prev
            else:
                Vinit[self.BtA] = self.Bt_prev[0]
                Vinit[self.BtB] = self.Bt_prev[1]
                Vinit[self.BtC] = self.Bt_prev[2]

            Vinit[self.lambda_BtA] = 0.0001
            Vinit[self.lambda_BtB] = 0.0001
            Vinit[self.lambda_BtC] = 0.0001

            Vinit[self.nodeA_dual_ineq_Bt] = 0.0001
            Vinit[self.nodeA_dual_ineq_Bt] = 0.0001
            Vinit[self.nodeA_dual_ineq_Bt] = 0.0001

            Vinit[self.nodeA_dual_ineq_Bt_upper] = 0.0001
            Vinit[self.nodeA_dual_ineq_Bt_upper] = 0.0001
            Vinit[self.nodeA_dual_ineq_Bt_upper] = 0.0001
        else:
            Vinit[self.node_p_plus] = 0.0001

            Vinit[self.node_p_minus] = 500

            Vinit[self.node_dual_ineq_var_p_plus] = 0.0001

            Vinit[self.node_dual_ineq_var_p_plus_upper] = 0.0001

            Vinit[self.node_dual_ineq_var_p_minus] = 0.0001

            Vinit[self.node_dual_ineq_var_p_minus_upper] = 0.0001

            # SOC Initialization

            Vinit[self.Bt_nodes[0]] = self.Bt_prev
            Vinit[self.dual_Bt_nodes[0]] = 0.0001
            Vinit[self.mu_index_Bt[0]] = 0.0001
            Vinit[self.mu_index_Bt_upper[0]] = 0.0001
        
        return Vinit


    def assign_nodes(self, node, node_index):
        # So that batteries are isolated from Infeasibility currents and several can be placed at an individual node
        self.isTriplex = node.isTriplex

        if (self.isTriplex == False):
            if (self.single_phase == ""):
                # Three-Phase Battery

                self.nodeA_p_plus = node_index.__next__()
                self.nodeB_p_plus = node_index.__next__()
                self.nodeC_p_plus = node_index.__next__()

                self.nodeA_p_minus = node_index.__next__()
                self.nodeB_p_minus = node_index.__next__()
                self.nodeC_p_minus = node_index.__next__()

                self.nodeA_dual_ineq_var_p_plus = node_index.__next__()
                self.nodeB_dual_ineq_var_p_plus = node_index.__next__()
                self.nodeC_dual_ineq_var_p_plus = node_index.__next__()

                self.nodeA_dual_ineq_var_p_plus_upper = node_index.__next__()
                self.nodeB_dual_ineq_var_p_plus_upper = node_index.__next__()
                self.nodeC_dual_ineq_var_p_plus_upper = node_index.__next__()

                self.nodeA_dual_ineq_var_p_minus = node_index.__next__()
                self.nodeB_dual_ineq_var_p_minus = node_index.__next__()
                self.nodeC_dual_ineq_var_p_minus = node_index.__next__()

                self.nodeA_dual_ineq_var_p_minus_upper = node_index.__next__()
                self.nodeB_dual_ineq_var_p_minus_upper = node_index.__next__()
                self.nodeC_dual_ineq_var_p_minus_upper = node_index.__next__()


                # Power distribution Nodes +/-
                self.P_plus_nodes = [self.nodeA_p_plus, self.nodeB_p_plus,
                    self.nodeC_p_plus]
                self.P_minus_nodes = [self.nodeA_p_minus, self.nodeB_p_minus,
                        self.nodeC_p_minus]
                
                # Lower Inequality Nodes set
                self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_var_p_plus, self.nodeB_dual_ineq_var_p_plus,
                        self.nodeC_dual_ineq_var_p_plus]
                self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_var_p_minus, self.nodeB_dual_ineq_var_p_minus,
                        self.nodeC_dual_ineq_var_p_minus]
                
                # Upper Inequality Nodes set
                self.dual_ineq_r_plus_nodes_upper = [self.nodeA_dual_ineq_var_p_plus_upper, 
                                    self.nodeB_dual_ineq_var_p_plus_upper,self.nodeC_dual_ineq_var_p_plus_upper]
                self.dual_ineq_r_minus_nodes_upper = [self.nodeA_dual_ineq_var_p_minus_upper, 
                                    self.nodeB_dual_ineq_var_p_minus_upper,self.nodeC_dual_ineq_var_p_minus_upper]
                
                # SOC equality terms
                # SOC evolution equality constraints
                self.BtA = node_index.__next__()
                self.BtB = node_index.__next__()
                self.BtC = node_index.__next__()
                self.lambda_BtA = node_index.__next__()
                self.lambda_BtB = node_index.__next__()
                self.lambda_BtC = node_index.__next__()

                # Upper and Lower Inequality bounds for SOC
                # lowerer inequality bounds
                self.nodeA_dual_ineq_Bt = node_index.__next__()
                self.nodeB_dual_ineq_Bt = node_index.__next__()
                self.nodeC_dual_ineq_Bt = node_index.__next__()

                # upper inequality bounds
                self.nodeA_dual_ineq_Bt_upper = node_index.__next__()
                self.nodeB_dual_ineq_Bt_upper = node_index.__next__()
                self.nodeC_dual_ineq_Bt_upper = node_index.__next__()

                self.Bt_nodes = [self.BtA, self.BtB, self.BtC]
                self.dual_Bt_nodes = [self.lambda_BtA, self.lambda_BtB, self.lambda_BtC]

                # SOC dynamics full set of inequalities
                self.mu_index_Bt = [self.nodeA_dual_ineq_Bt, self.nodeB_dual_ineq_Bt, self.nodeC_dual_ineq_Bt]
                self.mu_index_Bt_upper = [self.nodeA_dual_ineq_Bt_upper, self.nodeB_dual_ineq_Bt_upper, self.nodeC_dual_ineq_Bt_upper]

            else:
                # Single-Phase Battery

                self.node_p_plus = node_index.__next__()

                self.node_p_minus = node_index.__next__()

                self.node_dual_ineq_var_p_plus = node_index.__next__()

                self.node_dual_ineq_var_p_plus_upper = node_index.__next__()

                self.node_dual_ineq_var_p_minus = node_index.__next__()

                self.node_dual_ineq_var_p_minus_upper = node_index.__next__()


                # Power distribution Nodes +/-
                self.P_plus_nodes = [self.node_p_plus]
                self.P_minus_nodes = [self.node_p_minus]
                
                # Lower Inequality Nodes set
                self.dual_ineq_r_plus_nodes = [self.node_dual_ineq_var_p_plus]
                self.dual_ineq_r_minus_nodes = [self.node_dual_ineq_var_p_minus]
                
                # Upper Inequality Nodes set
                self.dual_ineq_r_plus_nodes_upper = [self.node_dual_ineq_var_p_plus_upper]
                self.dual_ineq_r_minus_nodes_upper = [self.node_dual_ineq_var_p_minus_upper]


                # SOC equality terms
                # SOC evolution equality constraints
                self.Bt = node_index.__next__()
                self.lambda_Bt = node_index.__next__()

                # Upper and Lower Inequality bounds for SOC
                # lowerer inequality bounds
                self.node_dual_ineq_Bt = node_index.__next__()

                # upper inequality bounds
                self.node_dual_ineq_Bt_upper = node_index.__next__()

                self.Bt_nodes = [self.Bt]
                self.dual_Bt_nodes = [self.lambda_Bt]

                # SOC dynamics full set of inequalities
                self.mu_index_Bt = [self.node_dual_ineq_Bt]
                self.mu_index_Bt_upper = [self.node_dual_ineq_Bt_upper]


            # Power Distribution full set
            self.infeas_real_var_index = self.P_plus_nodes + self.P_minus_nodes
            self.mu_index = self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes
            self.mu_index_upper = self.dual_ineq_r_plus_nodes_upper + self.dual_ineq_r_minus_nodes_upper


            """
            self.P_plus_nodes = node.P_plus_nodes
            self.P_minus_nodes = node.P_minus_nodes
            # lower-bound nodes
            self.dual_ineq_r_plus_nodes = node.dual_ineq_r_plus_nodes
            self.dual_ineq_r_minus_nodes = node.dual_ineq_r_minus_nodes
            # upper-bound nodes
            self.dual_ineq_r_plus_nodes_upper = [node.nodeA_dual_ineq_var_p_plus_upper, node.nodeB_dual_ineq_var_p_plus_upper, node.nodeC_dual_ineq_var_p_plus_upper]
            self.dual_ineq_r_minus_nodes_upper = [node.nodeA_dual_ineq_var_p_minus_upper, node.nodeB_dual_ineq_var_p_minus_upper, node.nodeC_dual_ineq_var_p_minus_upper]
            """


            if (self.source_type == 'PQ'):
                # TODO: define Triplex node assignments
                """
                #self.Q_nodes = node.Q_nodes
                self.Q_plus_nodes = node.Q_plus_nodes
                self.Q_minus_nodes = node.Q_minus_nodes
                # lower-bound nodes
                self.dual_ineq_i_nodes = node.dual_ineq_i_nodes
                # upper-bound nodes
                # TODO: have to add these nodes and store them internally in battery object
                self.dual_ineq_i_nodes_upper = [node.nodeA_dual_ineq_var_p_plus_upper, node.nodeB_dual_ineq_var_p_plus_upper, node.nodeC_dual_ineq_var_p_plus_upper]
                """
                raise Exception("Battery placed at Triplex node [Functionality Not Implemented]")
        else:
            # TODO: define Triplex node assignments
            raise Exception("Battery placed at Triplex node [Functionality Not Implemented]")

        return node_index
    
    def show_power_output(self, V):
        print("Battery {}".format(self.ID))
        print("Charging power:")
        for i,pch in enumerate(self.P_plus_nodes):
            phase_list = ["A","B","C"]
            print(" Pch {}: {}".format(phase_list[i], V[pch]), end="")

        print("\nDischarging power:")
        for i,pd in enumerate(self.P_minus_nodes):
            phase_list = ["A","B","C"]
            print(" Pd {}: {}".format(phase_list[i], V[pd]), end="")
        print("")

        

            

    def calc_residuals(self, node, V, res_eqn, cs_tol):
        node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
        node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]

        node_Lr = node.dual_eq_var_r_nodes
        node_Li = node.dual_eq_var_i_nodes

        if (self.single_phase == ""):
            P_plus_nodes = self.P_plus_nodes
            P_minus_nodes = self.P_minus_nodes

            dual_ineq_r_plus_nodes = self.dual_ineq_r_plus_nodes
            dual_ineq_r_minus_nodes = self.dual_ineq_r_minus_nodes

            dual_ineq_r_plus_nodes_upper = self.dual_ineq_r_plus_nodes_upper
            dual_ineq_r_minus_nodes_upper = self.dual_ineq_r_minus_nodes_upper

            Bt_nodes = self.Bt_nodes
            dual_Bt_nodes = self.dual_Bt_nodes

            mu_index_Bt = self.mu_index_Bt
            mu_index_Bt_upper = self.mu_index_Bt_upper

            num_of_phases = len(P_plus_nodes)
            for index in range(num_of_phases):
                Vr = V[node_Vr[index]]
                Vi = V[node_Vi[index]]
                Lr = V[node_Lr[index]]
                Li = V[node_Li[index]]


                ### power residuals

                P_plus = V[P_plus_nodes[index]] 
                P_minus = V[P_minus_nodes[index]]
                Q = 0

                partials = self.calculate_partial_derivatives(Vr, Vi, (P_plus - P_minus), Q, calc_hessian = True)

                res_eqn[node_Vr[index]] += partials['Ir']
                res_eqn[node_Vi[index]] += partials['Ii']

                res_eqn[node_Lr[index]] += partials['dIr_dVr']*Lr + partials['dIi_dVr']*Li
                res_eqn[node_Li[index]] += partials['dIr_dVi']*Lr + partials['dIi_dVi']*Li

                ### power inequality constraint residuals

                # defining lower bound mu's
                mu_pch_lb = V[dual_ineq_r_plus_nodes[index]]
                mu_pd_lb = V[dual_ineq_r_minus_nodes[index]]

                # defining upper bound mu's
                mu_pch_ub = V[dual_ineq_r_plus_nodes_upper[index]]
                mu_pd_ub = V[dual_ineq_r_minus_nodes_upper[index]]

                # defining residuals for power 
                res_eqn[P_plus_nodes[index]] += partials['dIr_dP']*Lr + partials['dIi_dP']*Li + self.C_ch*P_plus - mu_pch_lb + mu_pch_ub
                res_eqn[P_minus_nodes[index]] += -partials['dIr_dP']*Lr - partials['dIi_dP']*Li + self.C_d - mu_pd_lb + mu_pd_ub

                # complementary slack residuals (lower bounds)
                res_eqn[dual_ineq_r_plus_nodes[index]] += -mu_pch_lb*P_plus + cs_tol
                res_eqn[dual_ineq_r_minus_nodes[index]] += -mu_pd_lb*P_minus + cs_tol

                # complementary slack residuals (upper bounds)
                res_eqn[dual_ineq_r_plus_nodes_upper[index]] += mu_pch_ub*(P_plus-self.P_max) + cs_tol
                res_eqn[dual_ineq_r_minus_nodes_upper[index]] += mu_pd_ub*(P_minus-self.P_max) + cs_tol

                ### battery soc equality constraint residuals

                Bt = V[Bt_nodes[index]]
                lambda_Bt = V[dual_Bt_nodes[index]]
                mu_low_Bt = V[mu_index_Bt[index]]
                mu_up_Bt =  V[mu_index_Bt_upper[index]]

                # battery soc constraint
                res_eqn[Bt_nodes[index]] += Bt - self.Bt_prev[index]  - self.Mch*P_plus + self.Md*P_minus

                # battery lambda constraint
                res_eqn[dual_Bt_nodes[index]] += lambda_Bt  - self.Lb_next[index] + mu_up_Bt - mu_low_Bt

                ### battery soc inequality constraint residuals

                # complementary slack residuals (lower bounds)
                res_eqn[mu_index_Bt[index]] += -mu_low_Bt*Bt + cs_tol

                # complementary slack residuals (upper bounds)
                res_eqn[mu_index_Bt_upper[index]] += mu_up_Bt*(Bt-self.SOC_max) + cs_tol
        else:
            print(colored(f'TODO: residuals for single-phase battery not implemented', 'red'))
            

        return res_eqn
