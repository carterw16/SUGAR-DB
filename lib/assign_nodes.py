from itertools import count


def assign_nodes(node_key, nodes, regulators, transformers, switches, fuses, reactors, ibdgs, settings, obj=2, source_type = None, stamp_slack_bus_flag = False, stamp_tplx_bus_sources = False, stamp_infeas_neutral_source = False, casetype = 0):
    # Heuristic flags
    voltage_unbalance_enabled = settings["voltage unbalance settings"]["enforce voltage unbalance"]
    # create the node index
    node_index_ = count(0)
    # assign nodes for all grid objects
    # if running an infeasibility analysis, 
    for ele in nodes:
        if not ele.isTriplex:
            node_index_ = ele.assign_nodes(node_index_, obj, source_type, stamp_slack_bus_flag, stamp_infeas_neutral_source, casetype)
        else:
            node_index_ = ele.assign_triplex_nodes(node_index_, obj, source_type, stamp_tplx_bus_sources, stamp_infeas_neutral_source)
        
        if not ele.isTriplex and voltage_unbalance_enabled:
            node_index_ = ele.assign_voltage_unbalance_nodes(node_index_)

    for ele in regulators:
        node_index_ = ele.assign_nodes(node_key, node_index_, nodes)
        ele.check_connections()

    for ele in transformers:
        if ele.connect_type not in [ele._SINGLE_PHASE, ele._CENTER_TAP]:
            node_index_ = ele.assign_nodes_3phase(node_key, node_index_, nodes)

        else:
            node_index_ = ele.assign_nodes_centertap(node_key, node_index_, nodes)

    for ele in switches:
        node_index_ = ele.assign_nodes(node_index_)

    for ele in fuses:
        node_index_ = ele.assign_nodes(node_index_)

    for ele in reactors:
        node_index_ = ele.assign_nodes(node_key, node_index_)

    for ele in ibdgs:
        node_index_ = ele.assign_nodes(node_index_)

    return node_index_, nodes, regulators, transformers, switches, fuses, reactors, ibdgs
