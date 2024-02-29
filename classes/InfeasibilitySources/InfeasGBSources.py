from types import MethodType
from classes.InfeasibilitySources.base import InfeasibilitySources
# from .InfeasibilityStamps import ()

class InfeasGBSources(InfeasibilitySources):
    def __init__(self, node, obj, obj_type, obj_scalar, source_type):
        InfeasibilitySources.__init__(self, node)
        self.source_type = source_type
        self.assign_nodes(node, obj)

        self.obj_scaling = obj_scalar

    def stampY(self, i, j, val, Ynlin_val, Ynlin_row, Ynlin_col, idx):
        if val == 0.:
            return idx
        Ynlin_val[idx] = val
        Ynlin_row[idx] = i
        Ynlin_col[idx] = j
        idx += 1
        return idx


    def stampJ(self, i, val, Jnlin_val, Jnlin_row, idx):
        if val == 0.:
            return idx
        Jnlin_val[idx] = val
        Jnlin_row[idx] = i
        idx += 1
        return idx
    
    def assign_nodes(self, node, obj):
        self.isTriplex = node.isTriplex

        if obj == 'L1':
            # TODO: Develop L1 capacity
            temp = 1
        elif obj == 'L2':
            if node.isTriplex == True:
                if self.source_type == 'GB':
                    self.node1_G = node.node1_G
                    self.node2_G = node.node2_G
                    #self.nodeN_G = node.nodeN_G
                self.node1_B = node.node1_B
                self.node2_B = node.node2_B
               # self.nodeN_B = node.nodeN_B
            else:
                if self.source_type == 'GB':
                    self.nodeA_G = node.nodeA_G
                    self.nodeB_G = node.nodeB_G
                    self.nodeC_G = node.nodeC_G
                    #self.nodeN_G = node.nodeN_G
                self.nodeA_B = node.nodeA_B
                self.nodeB_B = node.nodeB_B
                self.nodeC_B = node.nodeC_B
                #self.nodeN_B = node.nodeN_B

            self.B_nodes = node.B_nodes
            if self.source_type == 'GB':
                self.G_nodes = node.G_nodes
    
    def stamp_nonlinear_infeasibility_GB(self, infeas_source, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):   
        # Stamping the addition of the infeasibility source in the network equations in Vr, Vi nodes
        Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_equality_constraint(infeas_source, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
        # Stamping equations from (d-Lagrange)/(d-primal variable)
        Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J= self.stamp_stationarity_constraints(infeas_source, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
        
        return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J
    
    def stamp_equality_constraint(self, infeas_source, node, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J):
        if infeas_source.source_type == 'GB':
            G_nodes = infeas_source.G_nodes
        B_nodes = infeas_source.B_nodes
        
        if infeas_source.isTriplex:
            num_of_phases = 2

            node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
            node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
        else:
            num_of_phases = 3

            node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
            node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]

        for index in range(num_of_phases):
            B = V[B_nodes[index]]
            if infeas_source.source_type == 'GB':
                G = V[G_nodes[index]]

            Vr = V[node_Vr[index]]
            Vi = V[node_Vi[index]]

            ################## EQUALITY CONSTRAINTS #############################################
            # Real Equality Constraint => ... + Vi_across*Bslack = 0
            # Equality Constraint => ...  [(Vi)^k * Bslack + (Bslack)^k*Vi - (Vi)^k * (Bslack)^k]
            # From:
            if infeas_source.source_type == 'GB':
                idx_Y = self.stampY(node_Vr[index], G_nodes[index],  Vr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(node_Vr[index], node_Vr[index],  G, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_J = self.stampJ(node_Vr[index], -B*Vi + Vr*G, Jnlin_val, Jnlin_row, idx_J)
            else:
                idx_J = self.stampJ(node_Vr[index], -B*Vi, Jnlin_val, Jnlin_row, idx_J)
            idx_Y = self.stampY(node_Vr[index], B_nodes[index], -Vi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Vr[index], node_Vi[index], -B, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # Imaginary Inequality Constraint => ... - Vr*Bslack = 0
            # Imaginary Equality Constraint => ... - [(Vr)^k * Bslack + (Bslack)^k*Vr - (Vr)^k * (Bslack)^k]
            # From: 
            if infeas_source.source_type == 'GB':
                idx_Y = self.stampY(node_Vi[index], G_nodes[index], Vi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(node_Vi[index], node_Vi[index], G, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_J = self.stampJ(node_Vi[index], G*Vi + Vr*B, Jnlin_val, Jnlin_row, idx_J)
            else:
                idx_J = self.stampJ(node_Vi[index], Vr*B, Jnlin_val, Jnlin_row, idx_J)
            idx_Y = self.stampY(node_Vi[index], B_nodes[index], Vr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Vi[index], node_Vr[index], B, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
        
        return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J
    
    def stamp_stationarity_constraints(self, infeas_source, node, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J):
        if infeas_source.source_type == 'GB':
            G_nodes = infeas_source.G_nodes
        B_nodes = infeas_source.B_nodes

        node_Lr = node.dual_eq_var_r_nodes
        node_Li = node.dual_eq_var_i_nodes
        
        if infeas_source.isTriplex:
            num_of_phases = 2

            node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
            node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]

        else:
            num_of_phases = 3

            node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
            node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
        
        for index in range(num_of_phases):
            B = V[B_nodes[index]]
            if infeas_source.source_type == 'GB':
                G = V[G_nodes[index]]
            else: 
                G = 0

            Vr = V[node_Vr[index]]
            Vi = V[node_Vi[index]]

            Lr = V[node_Lr[index]]
            Li = V[node_Li[index]]

            partials = self.calculate_partial_derivatives(Vr, Vi, G, B, Lr, Li, calc_hessian = True)
            ###################### STATIONARITY CONSTRAINTS ############################################
            # Equations for primal variable B
            # dL/dB = cost_val*2*B*(Vr^2 + Vi^2)^2 - Lr*Vi + Li*Vr
            idx_Y = self.stampY(B_nodes[index], B_nodes[index], partials['d2L_dB2'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(B_nodes[index], node_Vr[index], partials['d2L_dBdVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(B_nodes[index], node_Vi[index], partials['d2L_dBdVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(B_nodes[index], node_Lr[index], partials['d2L_dBdLr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(B_nodes[index], node_Li[index], partials['d2L_dBdLi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_J = self.stampJ(B_nodes[index], partials['dL_dB_hist'], Jnlin_val, Jnlin_row, idx_J)

            # Equations for primal variable G
            if infeas_source.source_type == 'GB':
                # dL/dG = alpha*2*G*(Vr^2 + Vi^2)^2 + Lr*Vr + Li*Vi
                idx_Y = self.stampY(G_nodes[index], G_nodes[index], partials['d2L_dG2'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(G_nodes[index], node_Vr[index], partials['d2L_dGdVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(G_nodes[index], node_Vi[index], partials['d2L_dGdVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(G_nodes[index], node_Lr[index], partials['d2L_dGdLr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stampY(G_nodes[index], node_Li[index], partials['d2L_dGdLi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)		
                idx_J = self.stampJ(G_nodes[index], partials['dL_dG_hist'], Jnlin_val, Jnlin_row, idx_J)
            
            # Equations for primal variable Vr
            # dL/dVr_from [stamped to Lr node rows]
            idx_Y = self.stampY(node_Lr[index], B_nodes[index], partials['d2L_dVrdB'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            if infeas_source.source_type == 'GB':
                idx_Y = self.stampY(node_Lr[index], G_nodes[index], partials['d2L_dVrdG'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Lr[index], node_Vr[index], partials['d2L_dVr2'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Lr[index], node_Vi[index], partials['d2L_dVrdVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Lr[index], node_Lr[index], partials['d2L_dVrdLr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Lr[index], node_Li[index], partials['d2L_dVrdLi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_J = self.stampJ(node_Lr[index], partials['dL_dVr_hist'], Jnlin_val, Jnlin_row, idx_J)

            # Equations for primal variable Vi
            # dL/dVi_from [stamped to Li node rows]
            idx_Y = self.stampY(node_Li[index], B_nodes[index], partials['d2L_dVidB'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            if infeas_source.source_type == 'GB':
                idx_Y = self.stampY(node_Li[index], G_nodes[index], partials['d2L_dVidG'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Li[index], node_Vr[index], partials['d2L_dVidVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Li[index], node_Vi[index], partials['d2L_dVi2'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Li[index], node_Lr[index], partials['d2L_dVidLr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_Y = self.stampY(node_Li[index], node_Li[index], partials['d2L_dVidLi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            idx_J = self.stampJ(node_Li[index], partials['dL_dVi_hist'], Jnlin_val, Jnlin_row, idx_J)
    
        return(Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J)
    
    def calculate_partial_derivatives(self, Vr, Vi, B, G, Lr, Li, calc_hessian):
        cost_val = self.obj_scaling
        partials = {}

        Vr2Vi2 = Vr**2 + Vi**2
        Vr2Vi2Sq = Vr2Vi2**2
        B2G2 = B**2 + G**2 

        ###################### PARTIAL DERIVATIVES ############################################
        # NOTE: written for L2 objective w/ alpha scaling parameter
        # dL/dB = alpha*2*B*(Vr^2 + Vi^2)^2 - Lr*Vi + Li*Vr
        partials['dL_dB'] = cost_val*2*B*Vr2Vi2Sq - Lr*Vi + Li*Vr
        partials['d2L_dB2'] = cost_val*2*Vr2Vi2Sq
        partials['d2L_dBdVr'] = cost_val*8*B*Vr2Vi2*Vr + Li
        partials['d2L_dBdVi'] = cost_val*8*B*Vr2Vi2*Vi - Lr
        partials['d2L_dBdLr'] = -Vi
        partials['d2L_dBdLi'] =  Vr
        partials['dL_dB_hist'] = partials['d2L_dB2']*B + partials['d2L_dBdVr']*Vr + partials['d2L_dBdVi']*Vi + partials['d2L_dBdLr']*Lr + partials['d2L_dBdLi']*Li - partials['dL_dB']

        # dL/dG = alpha*2*G*(Vr^2 + Vi^2)^2 + Lr*Vr + Li*Vi
        partials['dL_dG'] = cost_val*2*G*Vr2Vi2Sq + Lr*Vr + Li*Vi
        partials['d2L_dG2'] = cost_val*2*Vr2Vi2Sq
        partials['d2L_dGdVr'] = cost_val*8*G*Vr2Vi2*Vr + Lr
        partials['d2L_dGdVi'] = cost_val*8*G*Vr2Vi2*Vi + Li
        partials['d2L_dGdLr'] = Vr
        partials['d2L_dGdLi'] = Vi
        partials['dL_dG_hist'] = partials['d2L_dG2']*G + partials['d2L_dGdVr']*Vr + partials['d2L_dGdVi']*Vi + partials['d2L_dGdLr']*Lr + partials['d2L_dGdLi']*Li - partials['dL_dG']

        # dL/dVr = alpha*2*Vr*(G^2 + B^2)*(Vr^2 + Vi^2) + Lr*G + Li*B
        partials['dL_dVr'] = cost_val*4*Vr*B2G2*Vr2Vi2 + Lr*G + Li*B
        partials['d2L_dVrdB'] = cost_val*8*B*Vr*Vr2Vi2 + Li
        partials['d2L_dVrdG'] = cost_val*8*G*Vr*Vr2Vi2 + Lr
        partials['d2L_dVr2'] = cost_val*4*B2G2*(3*Vr**2 + Vi**2)
        partials['d2L_dVrdVi'] = cost_val*8*Vi*Vr*B2G2
        partials['d2L_dVrdLr'] = G
        partials['d2L_dVrdLi'] = B
        partials['dL_dVr_hist'] = partials['d2L_dVrdB']*B + partials['d2L_dVrdG']*G + partials['d2L_dVr2']*Vr + partials['d2L_dVrdVi']*Vi + partials['d2L_dVrdLr']*Lr + partials['d2L_dVrdLi']*Li - partials['dL_dVr']

        # dL/dVi = alpha*2*Vi*(G^2 + B^2)*(Vr^2 + Vi^2) - Lr*B + Li*G
        partials['dL_dVi'] = cost_val*4*Vi*B2G2*Vr2Vi2 - Lr*B + Li*G
        partials['d2L_dVidB'] = cost_val*8*B*Vi*Vr2Vi2 - Lr
        partials['d2L_dVidG'] = cost_val*8*G*Vr*Vr2Vi2 + Li
        partials['d2L_dVi2'] = cost_val*4*B2G2*(3*Vi**2 + Vr**2)
        partials['d2L_dVidVr'] = cost_val*8*B2G2*Vr*Vi
        partials['d2L_dVidLr'] = -B
        partials['d2L_dVidLi'] = G
        partials['dL_dVi_hist'] = partials['d2L_dVidB']*B + partials['d2L_dVidG']*G + partials['d2L_dVidVr']*Vr + partials['d2L_dVi2']*Vi + partials['d2L_dVidLr']*Lr + partials['d2L_dVidLi']*Li - partials['dL_dVi']

        return partials

