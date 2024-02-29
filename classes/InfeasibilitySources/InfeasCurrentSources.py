from types import MethodType
from classes.GlobalVars import _PHASE
from classes.InfeasibilitySources.base import InfeasibilitySources
from .InfeasibilityStamps import (stamp_ineq_duals,
                                  stamp_infeas_current_L2,
                                  stamp_infeas_current_L1, 
                                  stamp_infeas_dual_relationship, 
                                  stamp_linear_infeasibility, stamp_nonlinear_infeasibility,
                                  stamp_equality_constraint, stamp_stationarity_constraints)


class InfeasCurrentSources(InfeasibilitySources):
    def __init__(self, node, obj, obj_type, obj_scalar, source_type, stamp_neutral_infeas):
        InfeasibilitySources.__init__(self, node)  
        self.obj = obj
        self.obj_type = obj_type
        self.obj_scaling = obj_scalar
        self.source_type = source_type
        self.stamp_infeas_dual_relationship = MethodType(stamp_infeas_dual_relationship, self)
        self.stamp_linear_infeasibility = MethodType(stamp_linear_infeasibility, self)
        self.stamp_nonlinear_infeasibility = MethodType(stamp_nonlinear_infeasibility, self)
        if obj == 'L1' and self.obj_type == 'current':
            # # Infeasibility Stamps - Link Cython Methods to Class # #
            self.stamp_ineq_duals = MethodType(stamp_ineq_duals, self) 
            self.stamp_infeas_current_L1 = MethodType(stamp_infeas_current_L1, self)
       
        elif obj == 'L2' and self.obj_type == 'current':
            # # Infeasibility Stamps - Link Cython Methods to Class # #
            self.stamp_infeas_current_L2 = MethodType(stamp_infeas_current_L2, self)
        
        elif self.obj_type == 'power':
            self.stamp_equality_constraint = MethodType(stamp_equality_constraint, self)
            self.stamp_stationarity_constraints = MethodType(stamp_stationarity_constraints, self)
        
        self.assign_nodes(node, obj, stamp_neutral_infeas)

    def assign_nodes(self, node, obj, stamp_neutral_infeas):
        self.isTriplex = node.isTriplex

        if obj == 'L1':
            if node.isTriplex == True:
                self.node1_if_r_plus = node.node1_if_r_plus
                self.node1_if_r_minus = node.node1_if_r_minus
                self.node1_if_i_plus = node.node1_if_i_plus
                self.node1_if_i_minus = node.node1_if_i_minus

                self.node2_if_r_plus  = node.node2_if_r_plus 
                self.node2_if_r_minus = node.node2_if_r_minus
                self.node2_if_i_plus = node.node2_if_i_plus
                self.node2_if_i_minus = node.node2_if_i_minus

                self.node1_dual_ineq_r_minus = node.node1_dual_ineq_r_minus
                self.node1_dual_ineq_i_minus = node.node1_dual_ineq_i_minus
                self.node2_dual_ineq_r_minus = node.node2_dual_ineq_r_minus
                self.node2_dual_ineq_i_minus = node.node2_dual_ineq_i_minus
                
                self.node1_dual_ineq_r_plus = node.node1_dual_ineq_r_plus
                self.node1_dual_ineq_i_plus = node.node1_dual_ineq_i_plus
                self.node2_dual_ineq_r_plus = node.node2_dual_ineq_r_plus
                self.node2_dual_ineq_i_plus = node.node2_dual_ineq_i_plus

                if stamp_neutral_infeas:
                    self.nodeN_if_r_plus = node.nodeN_if_r_plus
                    self.nodeN_if_r_minus = node.nodeN_if_r_minus
                    self.nodeN_if_i_plus = node.nodeN_if_i_plus
                    self.nodeN_if_i_minus = node.nodeN_if_i_minus

                    self.nodeN_dual_ineq_r_minus = node.nodeN_dual_ineq_r_minus
                    self.nodeN_dual_ineq_i_minus = node.nodeN_dual_ineq_i_minus

                    self.nodeN_dual_ineq_r_plus = node.nodeN_dual_ineq_r_plus
                    self.nodeN_dual_ineq_i_plus = node.nodeN_dual_ineq_i_plus                 

            else:
                self.nodeA_if_r_plus = node.nodeA_if_r_plus
                self.nodeA_if_r_minus = node.nodeA_if_r_minus
                self.nodeA_if_i_plus = node.nodeA_if_i_plus
                self.nodeA_if_i_minus = node.nodeA_if_i_minus
                
                self.nodeB_if_r_plus = node.nodeB_if_r_plus
                self.nodeB_if_r_minus = node.nodeB_if_r_minus       
                self.nodeB_if_i_plus = node.nodeB_if_i_plus
                self.nodeB_if_i_minus = node.nodeB_if_i_minus
                
                self.nodeC_if_r_plus = node.nodeC_if_r_plus
                self.nodeC_if_r_minus = node.nodeC_if_r_minus            
                self.nodeC_if_i_plus = node.nodeC_if_i_plus
                self.nodeC_if_i_minus = node.nodeC_if_i_minus

                self.nodeA_dual_ineq_r_minus = node.nodeA_dual_ineq_r_minus
                self.nodeA_dual_ineq_i_minus = node.nodeA_dual_ineq_i_minus
                self.nodeB_dual_ineq_r_minus = node.nodeB_dual_ineq_r_minus
                self.nodeB_dual_ineq_i_minus = node.nodeB_dual_ineq_i_minus
                self.nodeC_dual_ineq_r_minus = node.nodeC_dual_ineq_r_minus
                self.nodeC_dual_ineq_i_minus = node.nodeC_dual_ineq_i_minus

                self.nodeA_dual_ineq_r_plus = node.nodeA_dual_ineq_r_plus
                self.nodeA_dual_ineq_i_plus = node.nodeA_dual_ineq_i_plus
                self.nodeB_dual_ineq_r_plus = node.nodeB_dual_ineq_r_plus
                self.nodeB_dual_ineq_i_plus = node.nodeB_dual_ineq_i_plus
                self.nodeC_dual_ineq_r_plus = node.nodeC_dual_ineq_r_plus
                self.nodeC_dual_ineq_i_plus = node.nodeC_dual_ineq_i_plus

                if stamp_neutral_infeas:
                    self.nodeN_if_r_plus = node.nodeN_if_r_plus
                    self.nodeN_if_r_minus = node.nodeN_if_r_minus             
                    self.nodeN_if_i_plus = node.nodeN_if_i_plus
                    self.nodeN_if_i_minus = node.nodeN_if_i_minus

                    self.nodeN_dual_ineq_r_minus = node.nodeN_dual_ineq_r_minus
                    self.nodeN_dual_ineq_i_minus = node.nodeN_dual_ineq_i_minus
                
                    self.nodeN_dual_ineq_r_plus = node.nodeN_dual_ineq_r_plus
                    self.nodeN_dual_ineq_i_plus = node.nodeN_dual_ineq_i_plus 

            self.dual_ineq_r_minus_nodes = node.dual_ineq_r_minus_nodes
            self.dual_ineq_r_plus_nodes = node.dual_ineq_r_plus_nodes
            self.dual_ineq_i_minus_nodes = node.dual_ineq_i_minus_nodes
            self.dual_ineq_i_plus_nodes = node.dual_ineq_i_plus_nodes
            
            self.if_r_plus_nodes = node.if_r_plus_nodes
            self.if_r_minus_nodes = node.if_r_minus_nodes
            self.if_i_plus_nodes = node.if_i_plus_nodes
            self.if_i_minus_nodes = node.if_i_minus_nodes

        elif obj == 'L2':
            if node.isTriplex == True:
                self.node1_if_r = node.node1_if_r
                self.node1_if_i = node.node1_if_i
                self.node2_if_r = node.node2_if_r
                self.node2_if_i = node.node2_if_i
                if stamp_neutral_infeas:
                    self.nodeN_if_r = node.nodeN_if_r
                    self.nodeN_if_i = node.nodeN_if_i
            else:
                self.nodeA_if_r = node.nodeA_if_r
                self.nodeA_if_i = node.nodeA_if_i
                self.nodeB_if_r = node.nodeB_if_r
                self.nodeB_if_i = node.nodeB_if_i
                self.nodeC_if_r = node.nodeC_if_r
                self.nodeC_if_i = node.nodeC_if_i
                if stamp_neutral_infeas:
                    self.nodeN_if_r = node.nodeN_if_r
                    self.nodeN_if_i = node.nodeN_if_i

            self.if_r_nodes = node.if_r_nodes
            self.if_i_nodes = node.if_i_nodes

    def calc_residuals_L1(self, V_out, node, epsilon, res_eqn):
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
        prim_flag = False
        dual_flag = False

        num_of_phases = len(self.dual_ineq_r_minus_nodes)
                        
        for i in range(0, num_of_phases):
            # Complementary Slackness
            res_eqn[self.dual_ineq_r_minus_nodes[i]] += V_out[self.dual_ineq_r_minus_nodes[i]] * V_out[node.if_r_minus_nodes[i]] - epsilon
            res_eqn[self.dual_ineq_i_minus_nodes[i]] += V_out[self.dual_ineq_i_minus_nodes[i]] * V_out[node.if_i_minus_nodes[i]] - epsilon           
            res_eqn[self.dual_ineq_r_plus_nodes[i]] += V_out[self.dual_ineq_r_plus_nodes[i]] * V_out[node.if_r_plus_nodes[i]] - epsilon
            res_eqn[self.dual_ineq_i_plus_nodes[i]] += V_out[self.dual_ineq_i_plus_nodes[i]] * V_out[node.if_i_plus_nodes[i]] - epsilon

            if prim_flag == False: 
                if (V_out[node.if_r_minus_nodes[i]] >= 0 and V_out[node.if_i_minus_nodes[i]] >= 0 
                    and V_out[node.if_r_plus_nodes[i]] >= 0 and V_out[node.if_i_plus_nodes[i]] >= 0):
                    prim_flag = False
                else:
                    prim_flag = True      
            
            if dual_flag == False:
                if (V_out[self.dual_ineq_r_minus_nodes[i]] >= 0 and V_out[self.dual_ineq_i_minus_nodes[i]] >= 0
                    and V_out[self.dual_ineq_r_plus_nodes[i]] >= 0 and V_out[self.dual_ineq_i_plus_nodes[i]] >= 0):
                    dual_flag = False
                else:
                    dual_flag = True
                                                    

        # L1 and power
        #elif self.infeas_settings['obj'] == 'L1 - power':
        # L2 and current
        #elif self.infeas_settings['obj'] == 'L2 - current':
        # L2 and power
        # elif self.infeas_settings['obj'] == 'L2 - power':

        return res_eqn, prim_flag, dual_flag