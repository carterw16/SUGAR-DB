"""Parses a dictionary with the power grid elements.

Author(s): Naeem Turner-Bandele, Elizabeth Foster
Created Date: 06-30-2020
Updated Date: 6-23-2023
Email: nturnerb@cmu.edu
Status: Development

Parses the power grid elements, initializes their respective classes, and then returns a list for each element.

"""
import time
import numpy as np
import re
from types import SimpleNamespace
from itertools import count

from classes.Capacitors import Capacitors
from classes.Fuses import Fuses
from classes.Generators import Slack
from classes.GlobalVars import object_parsed
from classes.IBDGs import FPNSC
from classes.InfeasibilitySources.InfeasCurrentSources import InfeasCurrentSources
from classes.InfeasibilitySources.InfeasGBSources import InfeasGBSources
from classes.InfeasibilitySources.InfeasPQSources import InfeasPQSources
from classes.InfeasibilitySources.BatterySources import BatterySources
from classes.lines.overhead_line import OverheadLines
from classes.lines.underground_lines import UndergroundLines
from classes.lines.triplex_lines import TriplexLines
from classes.Loads import ZIP, Constant, Exponential
from classes.Nodes import Nodes
from classes.Reactors import Reactors
from classes.Regulators import FixedTapRegulator, VariableRegulator
from classes.Switches import Switches
from classes.Transformers import Transformer
from classes.TriplexLoads import TriplexLoads

# Import Reader and Attribute Retriever
from lib.reader import reader
from lib.retrive_attributes import retrive_attributes
from lib.identify_voltage_type import identify_voltage_type
from lib.assign_nodes import assign_nodes

import ipdb


def parser(casename, settings, features = None):
    """Parses a three-phase test case.

        Calls the reader to get the data from the specified test cases. Then, uses retrieve_attributes to parse the
        case data from the reader. After, retrieving the attributes the objects are initialized using the properties
        of their respective classes.

        Args:
          casename (string):
            The name of the three-phase test case.
          features (dict):
            Indicates what special features are enabled and disabled.


        Returns:
          Lists with all initialized power grid elements including nodes, slacks, loads, triplex loads (tplxload),
          transformers (transformer), overhead lines (ohline), underground lines (ugline), triplex lines (tplxline),
          capacitors, regulators, switches, fuses, reactors, and distributed generators.
        """

    stamp_dual = settings["Stamp Dual"]
    if stamp_dual:
        obj = settings["infeas settings"]["obj"]
        obj_type = settings["infeas settings"]["obj type"]
        obj_scalar = settings["infeas settings"]["obj scalar"]
        try:
            infeas_node_list = settings["infeas settings"]["infeas node list"] 
        except KeyError:
            infeas_node_list = None

        # feature: battery sources
        try:
            battery_node_list = settings["infeas settings"]["battery_node_list"]
        except KeyError:
            battery_node_list = None
        

        try:
            stamp_slack_bus = settings["infeas settings"]["stamp slack bus"]
        except KeyError:
            stamp_slack_bus = True
        
        try:
            stamp_tplx_bus_sources = settings["infeas settings"]["triplex sources"]
        except KeyError:
            stamp_tplx_bus_sources = True
        
        try:
            stamp_neutral_infeas_source_flag = settings["infeas settings"]["neutral infeas source"]
        except KeyError:
            stamp_neutral_infeas_source_flag = False
       
        source_type = settings["infeas settings"]["source type"]
    else:
        obj = None
        stamp_slack_bus = False
        source_type = None
        stamp_tplx_bus_sources = False
        stamp_neutral_infeas_source_flag = False

    # Run the reader and retrive the case data
    casedata, simulation_stats = reader(casename)

    # Start the parser timer
    start_parser_time = time.time()

    # Identify the base voltage type based on the case type
    casetype = casedata["type"]

    # create the node_key and node_id iterators
    node_key = {}
    node_id_ = count(0)

    # # Parse node data and slack data # #
    # - Node Class Properties - #
    # ID | name | phases | nominal_voltage | node_connect_type | isTriplex |  bustype |
    # maxvoltage_error | busflags | referencebus  | meanrepairtime | degree | stamp_dual | voltage_type
    # - Slack Class Properties - #
    # ID | name | phases | nominal_voltage | phase_angle | connect_type | stamp_dual | voltage_type
    slack = []
    nodes = []
    slack_data = casedata['slacks'] if 'slacks' in casedata else None
    node_data = casedata['nodes']
        
    if slack_data:
        for key, data in slack_data.items():
            (_id, phases, nom_voltage, phase_angle, bustype, is_triplex,
             connect_type, voltage_setpoint, voltage_A, voltage_B,
             voltage_C, slack_parameters) = retrive_attributes(data, 'SWING', casetype)

            slack.append(
                Slack(_id,
                      _id,
                      phases,
                      nom_voltage,
                      phase_angle,
                      voltage_setpoint,
                      connect_type,
                      slack_parameters,
                      voltage_A=voltage_A,
                      voltage_B=voltage_B,
                      voltage_C=voltage_C,
                      stamp_dual = stamp_dual))
            if _id not in node_key:
                voltage_type, nom_voltage = identify_voltage_type(nom_voltage)
                nodes.append(
                    Nodes(_id,
                          _id,
                          phases,
                          nom_voltage,
                          voltage_setpoint,
                          connect_type,
                          is_triplex,
                          bustype=bustype,
                          voltage_type=voltage_type,
                          stamp_dual = stamp_dual))
                node_key[_id] = node_id_.__next__()

    for key, data in node_data.items():
        object_parsed.add('node')

        _id, phases, nom_voltage, bustype, is_triplex, connect_type, voltage_setpoint = retrive_attributes(
            data, 'node', casetype)

        voltage_type, nom_voltage = identify_voltage_type(nom_voltage)
        
        if _id not in node_key:
            nodes.append(
                    Nodes(_id,
                          _id,
                          phases,
                          nom_voltage,
                          voltage_setpoint,
                          connect_type,
                          is_triplex,
                          bustype=bustype,
                          voltage_type=voltage_type,
                          stamp_dual = stamp_dual))
            node_key[_id] = node_id_.__next__()

    # # Parse PQ and ZIP loads # #
    # # - Load Properties - #
    #  ID | name | parent | phases | is_delta | load_class | enable_CM | constant_power_A |
    # constant_power_B | constant_power_C | constant_current_A | constant_current_B | constant_current_C |
    # constant_impedance_A | constant_impedance_B | constant_impedance_C | constant_powerAB | constant_powerBC |
    # constant_power_CA

    load_data = casedata["loads"]
    load = []
    if load_data:
        for key, data in load_data.items():
            (_id, phases, parent, nom_voltage, is_delta, load_class, powerA,
             powerB, powerC, impedanceA, impedanceB, impedanceC, currentA,
             currentB, currentC, load_type, percent_Z, percent_I, percent_P,
             powerAB, powerBC, powerCA, connect_type, CVRwatts,
             CVRvars) = retrive_attributes(data, 'load', casetype)

            # check if voltage is LL or LN and adjust nom voltage if needed
            voltage_type, nom_voltage = identify_voltage_type(nom_voltage)

            # Create new node if the connecting element is the load itself
            if key == parent:
                nodes.append(Nodes(
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          1.0,
                          connect_type,
                          data.is_center_tap,
                          bustype=1,
                          stamp_dual = stamp_dual,
                          voltage_type=voltage_type
                          ))
                node_key[_id] = node_id_.__next__()
            if casetype == 1 and (parent not in node_key):
                nodes.append(Nodes(parent, parent, phases, nom_voltage, 1.0, connect_type, data.is_center_tap, bustype = 1, stamp_dual = stamp_dual, voltage_type = voltage_type))
                node_key[parent] = node_id_.__next__()
            if load_type == "ZIP":
                load.append(
                    ZIP(_id, _id, parent, phases, is_delta,
                        load_class, enable_CM, powerA, powerB, powerC,
                        percent_Z, percent_I, percent_P))
            elif load_type == "Exponential":
                load.append(
                    Exponential(_id, _id, parent, phases, is_delta,
                                nom_voltage, powerA, powerB, powerC, powerAB,
                                powerBC, powerCA, CVRwatts, CVRvars))
            else:
                load.append(
                    Constant(_id, _id, parent, phases, is_delta,
                             load_class, powerA, powerB, powerC,
                             currentA, currentB, currentC, impedanceA,
                             impedanceB, impedanceC, powerAB, powerBC, powerCA, nom_voltage, stamp_dual = stamp_dual))

    # # Parse Triplex Loads # #
    # # - Triplex Load Properties - #
    # name | ID | parent | phases | current1 | current12 | current2 | currentN | impedance1 | impedance2 |
    # nominal_voltage | power1 | power2 | shunt1 | shunt12 | shunt2
    tplxload_data = casedata[
        "triplex_loads"] if 'triplex_loads' in casedata else None
    tplxload = []
    if tplxload_data:
        for key, data in tplxload_data.items():
            (_id, parent, phases, current1, current12, current2, currentN,
             impedance1, impedance2, impedance12, nom_voltage, power1, power12,
             power2, shunt1, shunt12, shunt2,
             connect_type) = retrive_attributes(data, 'triplex_load', casetype)

            # check if voltage is LL or LN and adjust nom voltage if needed
            voltage_type, nom_voltage = identify_voltage_type(nom_voltage)

            if key == parent:
                nodes.append(Nodes(
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          connect_type,
                          data.is_center_tap,
                          bustype=1,
                          stamp_dual = stamp_dual,
                          voltage_type=voltage_type
                          ))
                node_key[_id] = node_id_.__next__()

            tplxload.append(
                TriplexLoads(_id, _id, parent, phases, current1, current12,
                             current2, currentN, impedance1, impedance12,
                             impedance2, nom_voltage, power1, power12, power2,
                             shunt1, shunt12, shunt2, stamp_dual = stamp_dual))
       
        
    # # Parse lines # #
    # # - Line Properties - # #
    # ID | name | type | nominal_voltage | from_node | to_node | length | impedance_matrix |
    # capacitance_matrix | spacing | unused_phases | conductor_phases
    # # - Overhead Line Properties - # #
    # Line Properties + | diameter
    # # - Underground Line Properties - # #
    # Line Properties + | outer_diameter | neutral_diameter | neutral_strands | conductor_diameter |
    # conductor_resistance | conductor_gmr | neutral_resistance | neutral_gmr

    ohline = []
    ugline = []
    tplxline = []
    line_data = casedata['lines']
    if line_data:
        for key, data in line_data.items():

            (_id, line_type, phases, nom_voltage, from_node, to_node, length,
             spacing, diameter, impedance_matrix, capacitance_matrix,
             outer_diameter, neutral_diameter, neutral_strands,
             conductor_diameter, conductor_resistance, conductor_gmr,
             neutral_resistance, neutral_gmr, unused_phases,
             conductor_phases, insulation_thickness, line_parameters) = retrive_attributes(data, 'line', casetype)

            if line_type == 'underground':
                ugline.append(
                    UndergroundLines(_id, key, line_type, phases, nom_voltage,
                                     from_node, to_node, line_parameters, length,
                                     impedance_matrix, capacitance_matrix,
                                     spacing, unused_phases, conductor_phases,
                                     outer_diameter, neutral_diameter,
                                     neutral_strands, conductor_diameter,
                                     conductor_resistance, conductor_gmr,
                                     neutral_resistance, neutral_gmr, stamp_dual = stamp_dual))
            elif line_type == 'triplex':
                tplxline.append(
                    TriplexLines(_id, key, line_type, phases, nom_voltage,
                                 from_node, to_node, line_parameters, length, impedance_matrix,
                                 capacitance_matrix, diameter, conductor_resistance, conductor_gmr, insulation_thickness, stamp_dual = stamp_dual))
            else:
                ohline.append(
                    OverheadLines(_id, key, line_type, phases, nom_voltage,
                                  from_node, to_node, line_parameters, length, spacing, unused_phases,
                                  diameter, conductor_phases, conductor_gmr, 
                                  conductor_resistance, impedance_matrix, capacitance_matrix,
                                  stamp_dual = stamp_dual))                                    
                                                                
            
    # # Parse capacitors # #
    # - Capacitor properties - #
    # name | parent | ID | phases | nominal_voltage | connect_type | capacitorA | capacitorB | capacitorC | pt_phase
    # | delay | upper_limit | lower_limit
    capacitor = []
    capacitor_data = casedata['capacitors']
    if capacitor_data:
        for key, data in capacitor_data.items():
            (_id, parent, phases, nom_voltage, connection_type, capacitorA,
             capacitorB, capacitorC, pt_phase, delay, upper_limit,
             lower_limit) = retrive_attributes(data, 'capacitor', casetype)

            voltage_type, nom_voltage = identify_voltage_type(nom_voltage)

            capacitor.append(
                Capacitors(_id, parent, _id, phases, nom_voltage,
                           connection_type, capacitorA, capacitorB, capacitorC,
                           pt_phase, delay, upper_limit, lower_limit, stamp_dual = stamp_dual))

    # # Parse regulators # #
    # - Regulator properties - #
    # name | ID | phases | from_node | to_node | connect_type | reg_type | tapA | tapB | tapC |
    # tapA_change_count | tapB_change_count | tapC_change_count | pt_phase | ct_phase
    regulator = []
    regulator_data = casedata['regulators']
    if regulator_data:
        for key, data in regulator_data.items():
            tap_controls = features["Tap Controls"] if "Tap Controls" in features else None

            (_id, phases, from_node, to_node, connection_type, reg_type, tapA,
             tapB, tapC, pt_phase, ct_phase, connected_transformer,
             connected_winding, bustype, voltage_setpoint, nom_voltage,
             band_center, bandwidth, pt_ratio, ct_prim, highstep, lowstep,
             rated_power, compensator_resistance, compensator_reactance,
             tap_status, aR_max, aR_min, num_tap_steps) = retrive_attributes(
                data, 'regulator', casetype, tap_controls)

            # check if voltage is LL or LN and adjust nom voltage if needed
            voltage_type, nom_voltage = identify_voltage_type(nom_voltage)

            # If transformer connected to regulator, create an intermediary connector node.
            if connected_transformer:
                if connected_winding == 1:
                    to_node = to_node + '_RegConnector'
                    flag_connector_exists = next(
                        (obj for obj in nodes if obj.ID == to_node), None)
                    if not flag_connector_exists:
                        nodes.append(Nodes(
                                  to_node,
                                  to_node,
                                  phases,
                                  nom_voltage,
                                  voltage_setpoint,
                                  connection_type,
                                  False,
                                  bustype=bustype,
                                  stamp_dual = stamp_dual,
                                  voltage_type=voltage_type))
                        node_key[_id] = node_id_.__next__()
                else:
                    from_node = from_node + '_RegConnector'
                    flag_connector_exists = next(
                        (obj for obj in nodes if obj.ID == from_node), None)
                    if not flag_connector_exists:
                        nodes.append(Nodes(
                                  from_node,
                                  from_node,
                                  phases,
                                  nom_voltage,
                                  voltage_setpoint,
                                  connection_type,
                                  False,
                                  bustype=bustype,
                                  stamp_dual = stamp_dual,
                                  voltage_type=voltage_type))
                        node_key[_id] = node_id_.__next__()

            # select fixed tap or variable regulators
            regulator_types = {"Fixed": FixedTapRegulator, "Variable": VariableRegulator}
            Regulator = regulator_types[tap_status]
            
            if tap_status == "Fixed":
                reg = Regulator(_id, _id, phases, from_node,
                              to_node, connection_type, reg_type, 
                              band_center, bandwidth, rated_power, 
                              lowstep, highstep, num_tap_steps, 
                              pt_ratio, pt_phase, ct_phase,
                              connected_transformer, connected_winding,
                              tapA, tapB, tapC, stamp_dual=stamp_dual)
            else: 
                reg = Regulator(_id, _id, phases, from_node,
                              to_node, connection_type, reg_type, 
                              band_center, bandwidth, rated_power, 
                              lowstep, highstep, num_tap_steps, 
                              pt_ratio, pt_phase, ct_phase,
                              connected_transformer, connected_winding,
                              stamp_dual=stamp_dual)
            regulator.append(reg)

    # # Parse transformers # #
    # # - Transformer Properties - #
    # name | ID | phases | from_node| to_node | nominal_voltage | connect_type | install_type |
    # power_rating | primary_voltage | secondary_voltage | resistance | reactance | shunt_impedance |
    # is_center_tap |
    # is_grounded | full_load_loss | no_load_loss | phase_shift
    transformer = []
    transformer_data = casedata['transformers']
    if transformer_data:
        for key, data in transformer_data.items():
            (_id, phases, from_node, to_node, nom_voltage, primary_voltage,
             secondary_voltage, connection_type, install_type, is_center_tap,
             is_grounded, reactance, resistance, shunt_impedance, loadloss,
             noload_loss, phase_shift,
             power_rating) = retrive_attributes(data, 'transformer', casetype)
            
            if nom_voltage is not None:
                # check if voltage is LL or LN and adjust nom voltage if needed
                voltage_type, nom_voltage = identify_voltage_type(nom_voltage)

            # if a regulator is connected to a transformer
            # search the list of regulators to find the regulator that the transformer is connected to
            connected_regulator = next(
                (reg for reg in regulator if reg.connected_transformer == key),
                None)

            if connected_regulator:
                # If transformer connected to regulator, reassign from or to nodes.
                connected_winding = getattr(connected_regulator,
                                            "connected_winding")
                if connected_winding == 1:
                    from_node = getattr(connected_regulator, "to_node")
                else:
                    to_node = getattr(connected_regulator, "from_node")

            transformer.append(
                Transformer(_id, _id, casetype, phases, from_node,
                            to_node, nom_voltage, connection_type, install_type,
                            power_rating, primary_voltage, secondary_voltage,
                            resistance, reactance, shunt_impedance,
                            connected_regulator, is_center_tap, is_grounded,
                            loadloss, noload_loss, phase_shift, stamp_dual=stamp_dual)
            )

    # # Switch # #
    # - Switch Properties - #
    # name | ID | from_node | to_node | phases | phase_A_state | phase_B_state | phase_C_state
    switch = []
    switch_data = casedata['switches'] if 'switches' in casedata else None

    if switch_data:
        for key, data in switch_data.items():
            (_id, from_node, to_node, phases, phase_A_state, phase_B_state,
             phase_C_state, switch_parameters) = retrive_attributes(data, 'switch', casetype)
            if casetype == 1:
                switch_parameters['file type'] = 'DSS'
            else:
                switch_parameters['file type'] = 'GLD'

            switch.append(
                    Switches(_id, _id, from_node, to_node, phases, switch_parameters, 
                            phase_A_state, phase_B_state, phase_C_state, stamp_dual = stamp_dual))

    # # Parse fuses # #
    # - Fuse properties - #
    # name | ID | from_node | to_node | phases | phase_A_state | phase_B_state | phase_C_state
    fuse = []
    fuse_data = casedata['fuses'] if 'fuses' in casedata else None
    if fuse_data:
        for key, data in fuse_data.items():
            (_id, from_node, to_node, phases, phase_A_state, phase_B_state,
             phase_C_state) = retrive_attributes(data, 'fuse')

            fuse.append(
                Fuses(_id, _id, from_node, to_node, phases,
                      phase_A_state, phase_B_state, phase_C_state, stamp_dual = stamp_dual))

    # # Parse reactors # #
    # - Reactor properties - #
    # node | ID | name | from_node | to_node | connection_type | phases | resistance | reactance
    reactor = []
    reactor_data = casedata['reactors'] if 'reactors' in casedata else None
    if reactor_data:
        for key, data in reactor_data.items():
            _id, from_node, to_node, phases, connection_type, resistance, reactance = retrive_attributes(
                data, 'reactor')
            reactor.append(
                Reactors(nodes, _id, key, from_node, to_node, connection_type,
                         phases, resistance, reactance, stamp_dual = stamp_dual))

    # # Parse IBDGs # #
    # - IBDG properties - #
    ibdg = []
    ibdg_data = casedata[
        'photovoltaics'] if 'photovoltaics' in casedata else None
    connected_loads = []

    if ibdg_data:
        for key, data in ibdg_data.items():
            (_id, phases, nom_voltage, connection_type, rated_power,
             active_rating, reactive_rating, power_factor, pv_power_max,
             irradiance, is_triplex) = retrive_attributes(data, 'ibdg')

            # check if voltage is LL or LN and adjust nom voltage if needed
            nom_voltage = nom_voltage if nom_voltage is not None else nodes[node_key[_id]].Vnom
            voltage_type, nom_voltage = identify_voltage_type(nom_voltage)
            voltage_setpoint = nodes[node_key[_id]].Vset
            # # Find the Connected Load and Its Powers # #
            find_end = re.match('.+([0-9])[^0-9]*$', _id)
            last_idx = find_end.start(1)
            load_name = _id[:last_idx + 1]
            if 'TMP0009' in casename:
                load_id = next((ele.from_node for ele in transformer if load_name in ele.to_node), None)
                connected_load = next((ele for ele in load if ele.ID == load_id), None)
            else:
                connected_load = next((ele for ele in load if ele.ID == load_name), None)

            P = active_rating
            Q = reactive_rating
            if connected_load:
                connected_loads.append(connected_load.ID)
                if phases == 7 or phases == 15:
                    P = np.sum(np.asarray(connected_load._cP))
                    Q = np.sum(np.asarray(connected_load._cQ))
                elif phases == 1:
                    P = float(connected_load._cP[0])
                    Q = float(connected_load._cQ[0])
                elif phases == 2:
                    P = float(connected_load._cP[1])
                    Q = float(connected_load._cQ[1])
                elif phases == 4:
                    P = float(connected_load._cP[2])
                    Q = float(connected_load._cQ[2])

            P = np.around(P, 0)
            Q = np.around(Q, 0)

            ibdg_settings = features["IBDGs"] if "IBDGs" in features else None
            inverter_settings = ibdg_settings[
                "Settings"] if "Settings" in ibdg_settings else None

            ibdg.append(
                FPNSC(_id,
                      data.name,
                      phases,
                      nom_voltage,
                      voltage_setpoint,
                      connection_type,
                      rated_power,
                      P,
                      Q,
                      active_rating,
                      reactive_rating,
                      power_factor,
                      irradiance,
                      pv_power_max,
                      inverter_settings,
                      is_triplex,
                      fixed_Q=False,
                      sensor="GVGC"))
    
    # ASSIGN NODES
    # If doing an infeasibility analysis, the infeasibility node indices are assigned
    # by source type in this function in addition to the regular indice assignments
    # TODO: use infeas_node_list to add infeasibility source only at specific nodes
    (node_index_, nodes, regulator, transformer, switch,
         fuse, reactor, ibdg) = assign_nodes(node_key, nodes, regulator,
                                                                        transformer, switch,
                                                                        fuse, reactor,
                                                                        ibdg, settings, obj, 
                                                                        source_type, stamp_slack_bus,
                                                                        stamp_tplx_bus_sources,stamp_neutral_infeas_source_flag,
                                                                        casetype)

    infeasibility_sources = []
    if stamp_dual:
        if infeas_node_list:
            # TODO: add assigning only to specific node list
            temp_todo_ = 1
        else:
            print("Assigning infeasibility sources to all nodes in the system")
            if source_type == 'current':
                for index in range(len(nodes)):
                    if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                        if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                            continue
                        else:
                            infeasibility_sources.append(InfeasCurrentSources(nodes[index], obj, obj_type, obj_scalar, source_type, stamp_neutral_infeas_source_flag))
            elif source_type == 'PQ' or source_type == 'Q':
                for index in range(len(nodes)):
                    if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                        if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                            continue
                        else:
                            infeasibility_sources.append(InfeasPQSources(nodes[index], obj, obj_type, obj_scalar, source_type, stamp_neutral_infeas_source_flag))
            elif source_type == 'GB' or source_type == 'B':
                for index in range(len(nodes)):
                    if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                        if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                            continue
                        else:
                            infeasibility_sources.append(InfeasGBSources(nodes[index], obj, obj_type, obj_scalar, source_type, stamp_neutral_infeas_source_flag))


    battery_sources = []
    # Battery source settings obj_type = "P" or "PQ", each battery has similar objective
    if stamp_dual:
        if battery_node_list:
            for bat in battery_node_list:
                for index in range(len(nodes)):
                    if (bat["node"] == nodes[index].ID):
                        if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                            # TODO: may need to put batteries at triplex nodes
                            if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                                continue
                            else:
                                _obj_type = bat["type"]
                                bat = BatterySources(nodes[index], node_index_, obj_scalar, _obj_type, stamp_neutral_infeas_source_flag,
                                                    ID = bat["ID"], P_max = bat["P_max"], P_min = bat["P_min"],
                                                    Mch = bat["Mch"], Md = bat["Md"], Bt_prev = bat["Bt_prev"],
                                                    C_ch = bat["C_ch"], C_d = bat["C_d"], single_phase = bat["single_phase"])
                                # append Battery to battery sources list and update node index
                                battery_sources.append(bat)
                                node_index_ = bat.updated_node_index

        """                      
        else:
            print("!! Assigning Batteries to every node!!")
            if source_type == 'P':
                for index in range(len(nodes)):
                    if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                        # TODO: may need to put batteries at triplex nodes
                        if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                            continue
                        else:
                            infeasibility_sources.append(BatterySources(nodes[index], _obj_type, obj_scalar, source_type, stamp_neutral_infeas_source_flag))
            elif source_type == 'PQ':
                for index in range(len(nodes)):
                    if (nodes[index].bustype == 3 and stamp_slack_bus) or (nodes[index].bustype != 3):
                        # TODO: may need to put batteries at triplex nodes
                        if stamp_tplx_bus_sources == False and nodes[index].isTriplex == True:
                            continue
                        else:
                            infeasibility_sources.append(BatterySources(nodes[index], _obj_type, obj_scalar, source_type, stamp_neutral_infeas_source_flag))
        """
    
    end_parser_time = time.time()

    # compute total parsing time and store in simulation stats
    parser_time = end_parser_time - start_parser_time
    simulation_stats.append(parser_time)
    casedata = {
        'node': nodes,
        'slack': slack,
        'load': load,
        'triplex_load': tplxload,
        'xfmr': transformer,
        'ohline': ohline,
        'ugline': ugline,
        'triplex_line': tplxline,
        'capacitor': capacitor,
        'regulator': regulator,
        'switch': switch,
        'fuse': fuse,
        'reactor': reactor,
        'ibdg': ibdg,
        'stats': simulation_stats,
        'infeasibility_sources': infeasibility_sources,
        'battery_sources': battery_sources,
        'voltage_bounds': None
    }

    
    casedata = SimpleNamespace(**casedata)
    return casedata, node_key, node_index_
