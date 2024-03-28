#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 15:01:04 2021

@author: emfoster
"""

import numpy as np
from classes.Nodes import Nodes

def calc_infeas_obj(V, node, xfmr, regulator, switch):
    feas_curr = []
    for i in range(0, len(node)):
        if node[i].isTriplex:
            continue
        else:
            feas_curr.append(V[node[i].nodeA_dual_eq_var_r])
            feas_curr.append(V[node[i].nodeA_dual_eq_var_i])
            feas_curr.append(V[node[i].nodeB_dual_eq_var_r])
            feas_curr.append(V[node[i].nodeB_dual_eq_var_i])
            feas_curr.append(V[node[i].nodeC_dual_eq_var_r])
            feas_curr.append(V[node[i].nodeC_dual_eq_var_i])
            feas_curr.append(V[node[i].nodeN_dual_eq_var_r])
            feas_curr.append(V[node[i].nodeN_dual_eq_var_i])
          
            if node[i].bustype == 3:
                feas_curr.append(V[node[i].slack_nodeA_Lr])
                feas_curr.append(V[node[i].slack_nodeA_Li])
                feas_curr.append(V[node[i].slack_nodeB_Lr])
                feas_curr.append(V[node[i].slack_nodeB_Li])
                feas_curr.append(V[node[i].slack_nodeC_Lr])
                feas_curr.append(V[node[i].slack_nodeC_Li])
    
    for i in range(0, len(xfmr)):
        feas_curr.append(V[xfmr[i].nodeA_Lrp])
        feas_curr.append(V[xfmr[i].nodeA_Lip])
        feas_curr.append(V[xfmr[i].nodeB_Lrp])
        feas_curr.append(V[xfmr[i].nodeB_Lip])
        feas_curr.append(V[xfmr[i].nodeC_Lrp])
        feas_curr.append(V[xfmr[i].nodeC_Lip])
        feas_curr.append(V[xfmr[i].nodeGnd_Lrp])
        feas_curr.append(V[xfmr[i].nodeGnd_Lip])
        feas_curr.append(V[xfmr[i].nodeGnd_Lrs])
        feas_curr.append(V[xfmr[i].nodeGnd_Lis])
        feas_curr.append(V[xfmr[i].nodeXtraA_Lrs])
        feas_curr.append(V[xfmr[i].nodeXtraA_Lis])        
        feas_curr.append(V[xfmr[i].nodeXtraB_Lrs])
        feas_curr.append(V[xfmr[i].nodeXtraB_Lis])
        feas_curr.append(V[xfmr[i].nodeXtraC_Lrs])
        feas_curr.append(V[xfmr[i].nodeXtraC_Lis])
    
    for i in range(0, len(regulator)):
        feas_curr.append(V[regulator[i].nodeLrA_rp])
        feas_curr.append(V[regulator[i].nodeLiA_ip])
        feas_curr.append(V[regulator[i].nodeLrB_rp])
        feas_curr.append(V[regulator[i].nodeLiB_ip])
        feas_curr.append(V[regulator[i].nodeLrC_rp])
        feas_curr.append(V[regulator[i].nodeLiC_ip])
        feas_curr.append(V[regulator[i].nodeLrGnd_rp])
        feas_curr.append(V[regulator[i].nodeLiGnd_ip])
        feas_curr.append(V[regulator[i].nodeLrGnd_rs])
        feas_curr.append(V[regulator[i].nodeLiGnd_is])
        feas_curr.append(V[regulator[i].nodeLrXtraA_rs])
        feas_curr.append(V[regulator[i].nodeLiXtraA_is])        
        feas_curr.append(V[regulator[i].nodeLrXtraB_rs])
        feas_curr.append(V[regulator[i].nodeLiXtraB_is])
        feas_curr.append(V[regulator[i].nodeLrXtraC_rs])
        feas_curr.append(V[regulator[i].nodeLiXtraC_is])
    
    for i in range(0, len(switch)):
        feas_curr.append(V[switch[i].node_switchA_Lr])
        feas_curr.append(V[switch[i].node_switchA_Li])
        feas_curr.append(V[switch[i].node_switchB_Lr])
        feas_curr.append(V[switch[i].node_switchB_Li])
        feas_curr.append(V[switch[i].node_switchC_Lr])
        feas_curr.append(V[switch[i].node_switchC_Li])
        
    obj = 1/2*np.square(np.linalg.norm(feas_curr, 2))
    return obj