import numpy as np


def check_error(err_nodes, V, Vsol, err_max_store, nodes, regulators, ibdgs):
    Verr = V[err_nodes]
    Vsol_err = Vsol[err_nodes]
    err = Verr - Vsol_err
    err_max = np.amax(abs(err))
    err_argmax = err_nodes[np.argmax(abs(err))]
    err_max_store.append(err_max)

    err_max_node = [(ele.ID, 'Var') for ele in nodes if hasattr(ele, 'nodeA_Vr') and ele.nodeA_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vbr') for ele in nodes if hasattr(ele, 'nodeB_Vr') and ele.nodeB_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vcr') for ele in nodes if hasattr(ele, 'nodeC_Vr') and ele.nodeC_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vai') for ele in nodes if hasattr(ele, 'nodeA_Vi') and ele.nodeA_Vi == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vbi') for ele in nodes if hasattr(ele, 'nodeB_Vi') and ele.nodeB_Vi == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vci') for ele in nodes if hasattr(ele, 'nodeC_Vi') and ele.nodeC_Vi == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'V1r') for ele in nodes if hasattr(ele, 'node1_Vr') and ele.node1_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'V2r') for ele in nodes if hasattr(ele, 'node2_Vr') and ele.node2_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vnr') for ele in nodes if hasattr(ele, 'nodeN_Vr') and ele.nodeN_Vr == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'V1i') for ele in nodes if hasattr(ele, 'node1_Vi') and ele.node1_Vi == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'V2i') for ele in nodes if hasattr(ele, 'node2_Vi') and ele.node2_Vi == err_argmax]

    if not err_max_node:
        err_max_node = [(ele.ID, 'Vni') for ele in nodes if hasattr(ele, 'nodeN_Vi') and ele.nodeN_Vi == err_argmax]

    if ibdgs:
        if not err_max_node:
            err_max_node = [(ele.ID, 'IBDG_Q') for ele in ibdgs if
                            hasattr(ele, 'node_Q3') and ele.node_Q3 == err_argmax]

    if regulators:
        if not err_max_node:
            err_max_node = [(ele.ID, 'aR_A') for ele in regulators if hasattr(ele, 'nodeA_aR') and ele.nodeA_aR ==
                            err_argmax]
        if not err_max_node:
            err_max_node = [(ele.ID, 'Va_diff') for ele in regulators if hasattr(ele, 'nodeA_V_diff') and
                            ele.nodeA_V_diff == err_argmax]
        if not err_max_node:
            err_max_node = [(ele.ID, 'aR_B') for ele in regulators if hasattr(ele, 'nodeB_aR') and ele.nodeB_aR ==
                            err_argmax]
        if not err_max_node:
            err_max_node = [(ele.ID, 'Vb_diff') for ele in regulators if hasattr(ele, 'nodeB_V_diff') and
                            ele.nodeB_V_diff == err_argmax]
        if not err_max_node:
            err_max_node = [(ele.ID, 'aR_C') for ele in regulators if hasattr(ele, 'nodeC_aR') and ele.nodeC_aR ==
                            err_argmax]
        if not err_max_node:
            err_max_node = [(ele.ID, 'Vc_diff') for ele in regulators if hasattr(ele, 'nodeC_V_diff') and
                            ele.nodeC_V_diff == err_argmax]

    return err_max, err_max_store, err_max_node
