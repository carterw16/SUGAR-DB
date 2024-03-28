"""Parses a dictionary with the power grid elements.

Author(s): Naeem Turner-Bandele
Created Date: 06-30-2020
Updated Date: 07-23-2020
Email: nturnerb@cmu.edu
Status: Development

Parses the power grid elements, initializes their respective classes, and then returns a list for each element.

"""
# Import Classes
from classes.Capacitors import Capacitors
from classes.Fuses import Fuses
from classes.Generators import Slack
from classes.GlobalVars import node_index_, object_parsed
from classes.Lines import OverheadLines, TriplexLines, UndergroundLines
from classes.Loads import ZIP, Constant, Exponential
from classes.Nodes import Nodes
from classes.Reactors import Reactors
from classes.Regulators import Regulators
from classes.Switches import Switches
from classes.Transformers import Transformer
from classes.TriplexLoads import TriplexLoads
from lib.reader import reader
# Import Reader and Attribute Retriever
from lib.retrive_attributes import retrive_attributes


def parser(casename,
           stampdual,
           features,
           freq_flag=False,
           dg_type=None,
           dgs_flag=False,
           gen_setpoints=None,
           enable_CM=False):
    """Parses a three-phase test case.

        Calls the reader to get the data from the specified test cases. Then, uses retrieve_attributes to parse the
        case data from the reader. After, retrieving the attributes the objects are initialized using the properties
        of their respective classes.

        Args:
          casename:
            The name of the three-phase test case. Type string.
          stamp_dual:
            Specifies if an optimal power flow is being performed. Currently an unused variable. Type bool.
          dg_type:
            Indicates the type of inverter-based distributed generator if DG's are used. Type string.
          dgs_flag:
            Specifies is DG's are enabled. Type bool.
          enable_CM:
            Indicates if current measurements are being taken. Type bool.

        Returns:
          Lists with all initialized power grid elements including nodes, slacks, loads, triplex loads (tplxload),
          transformers (xfmr), overhead lines (ohline), underground lines (ugline), triplex lines (tplxline),
          capacitors, regulators, switches, fuses, reactors, and distributed generators.
        """

    # Run the reader and retrive the case data
    casedata = reader(casename)
    # Identify the base voltage type based on the case type
    casetype = casedata["type"]
    if casetype != 0:
        voltage_type = 'LL'
    else:
        voltage_type = 'LN'

    # # Parse node data and slack data # #
    # - Node Class Properties - #
    # node_index_ | ID | name | phases | nominal_voltage | node_connect_type | isTriplex |  bustype |
    # maxvoltage_error | busflags | referencebus  | meanrepairtime | degree | stamp_dual | voltage_type
    # - Slack Class Properties - #
    # ID | name | phases | nominal_voltage | connect_type | stamp_dual | voltage_type
    slack = []
    node = []
    slack_data = casedata['slacks'] if 'slacks' in casedata else None
    node_data = casedata['nodes']

    if slack_data:
        for key, data in slack_data.items():
            _id, phases, nom_voltage, bustype, is_triplex, \
            connect_type, voltage_setpoint = retrive_attributes(data, 'SWING')
            slack.append(
                Slack(_id,
                      _id,
                      phases,
                      nom_voltage,
                      voltage_setpoint,
                      connect_type,
                      voltage_type=voltage_type))
            if _id not in Nodes.nodeKey:
                node.append(
                    Nodes(node_index_,
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          voltage_setpoint,
                          connect_type,
                          is_triplex,
                          bustype=bustype,
                          voltage_type=voltage_type,
                          stamp_dual = stampdual))

    for key, data in node_data.items():
        object_parsed.add('node')

        _id, phases, nom_voltage, bustype, is_triplex, connect_type, voltage_setpoint = retrive_attributes(
            data, 'node')

        if _id not in Nodes.nodeKey:
            try:
                node.append(
                    Nodes(node_index_,
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          voltage_setpoint,
                          connect_type,
                          is_triplex,
                          bustype=bustype,
                          voltage_type=voltage_type,
                          stamp_dual = stampdual))
            except TypeError:
                pass

    # # Parse PQ and ZIP loads # #
    # # - Load Properties - #
    # node_index_ | ID | name | parent | phases | is_delta | load_class | enable_CM | constant_power_A |
    # constant_power_B | constant_power_C | constant_current_A | constant_current_B | constant_current_C |
    # constant_impedance_A | constant_impedance_B | constant_impedance_C | constant_powerAB | constant_powerBC |
    # constant_power_CA

    load_data = casedata["loads"]
    load = []
    if load_data:
        for key, data in load_data.items():
            _id, phases, parent, nom_voltage, is_delta, load_class, \
            powerA, powerB, powerC, impedanceA, impedanceB, impedanceC, \
            currentA, currentB, currentC, load_type, percent_Z, percent_I, percent_P, powerAB, powerBC, powerCA, \
            connect_type, CVRwatts, CVRvars = retrive_attributes(data, 'load')

            # Create new node if the connecting element is the load itself
            if key == parent:
                node.append(
                    Nodes(node_index_,
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          1.0,
                          connect_type,
                          data.is_center_tap,
                          bustype=1,
                          stamp_dual = stampdual,
                          voltage_type=voltage_type
                          ))

            if load_type == "ZIP":
                load.append(
                    ZIP(node_index_, _id, _id, parent, phases, is_delta,
                        load_class, enable_CM, powerA, powerB, powerC,
                        percent_Z, percent_I, percent_P))
            elif load_type == "Exponential":
                load.append(
                    Exponential(node_index_, _id, _id, parent, phases, is_delta,
                                nom_voltage, powerA, powerB, powerC, powerAB,
                                powerBC, powerCA, CVRwatts, CVRvars))
            else:
                load.append(
                    Constant(node_index_, _id, _id, parent, phases, stampdual, is_delta,
                             load_class, enable_CM, powerA, powerB, powerC,
                             currentA, currentB, currentC, impedanceA,
                             impedanceB, impedanceC, powerAB, powerBC, powerCA))

    # # Parse Triplex Loads # #
    # # - Triplex Load Properties - #
    # name | ID | parent | phases | current1 | current12 | current2 | currentN | impedance1 | impedance2 |
    # nominal_voltage | power1 | power2 | shunt1 | shunt12 | shunt2
    tplxload_data = casedata[
        "triplex_loads"] if 'triplex_loads' in casedata else None
    tplxload = []

    if tplxload_data:
        for key, data in tplxload_data.items():
            _id, parent, phases, current1, current12, current2, currentN, impedance1, impedance2, impedance12, \
            nom_voltage, power1, power12, power2, shunt1, shunt12, shunt2, connect_type = retrive_attributes(data,
                                                                                                    'triplex_load')

            if key == parent:
                node.append(
                    Nodes(node_index_,
                          _id,
                          _id,
                          phases,
                          nom_voltage,
                          connect_type,
                          data.is_center_tap,
                          bustype=1,
                          voltage_type=voltage_type,
                          stamp_dual = stampdual))

            tplxload.append(
                TriplexLoads(_id, _id, parent, phases, current1, current12,
                             current2, currentN, impedance1, impedance12,
                             impedance2, nom_voltage, power1, power12, power2,
                             shunt1, shunt12, shunt2))

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

            _id, line_type, phases, nom_voltage, from_node, to_node, length, spacing, diameter, impedance_matrix, \
            capacitance_matrix, outer_diameter, neutral_diameter, neutral_strands, conductor_diameter, \
            conductor_resistance, conductor_gmr, neutral_resistance, neutral_gmr, unused_phases, conductor_phases = \
             retrive_attributes(data, 'line', casetype)

            if line_type == 'underground':
                ugline.append(
                    UndergroundLines(_id, key, line_type, phases, nom_voltage,
                                     from_node, to_node, length,
                                     impedance_matrix, capacitance_matrix,
                                     spacing, unused_phases, conductor_phases,
                                     outer_diameter, neutral_diameter,
                                     neutral_strands, conductor_diameter,
                                     conductor_resistance, conductor_gmr,
                                     neutral_resistance, neutral_gmr))
            elif line_type == 'triplex':
                tplxline.append(
                    TriplexLines(_id, key, line_type, phases, nom_voltage,
                                 from_node, to_node, length, impedance_matrix,
                                 capacitance_matrix))
            else:
                ohline.append(
                    OverheadLines(_id, key, line_type, phases, nom_voltage,
                                  from_node, to_node, length, impedance_matrix,
                                  capacitance_matrix, spacing, unused_phases,
                                  diameter, conductor_phases))

    # # Parse transformers # #
    # # - Transformer Properties - #
    # node_index_| node | name | ID | phases | from_node| to_node | nominal_voltage | connect_type | install_type |
    # power_rating | primary_voltage | secondary_voltage | resistance | reactance | shunt_impedance | is_center_tap |
    # is_grounded | full_load_loss | no_load_loss | phase_shift
    xfmr = []
    xfmr_data = casedata['transformers']

    if xfmr_data:
        for key, data in xfmr_data.items():
            if 'reg' not in key:
                _id, phases, from_node, to_node, nom_voltage, primary_voltage, secondary_voltage, connection_type, \
                install_type, is_center_tap, is_grounded, reactance, resistance, shunt_impedance, loadloss, \
                noload_loss, phase_shift, power_rating = retrive_attributes(data, 'xfmr', casetype)

                xfmr.append(
                    Transformer(node_index_, node, _id, _id, phases, from_node,
                                to_node, nom_voltage, connection_type,
                                install_type, power_rating, primary_voltage,
                                secondary_voltage, resistance, reactance,
                                shunt_impedance, is_center_tap, is_grounded,
                                loadloss, noload_loss, phase_shift, stamp_dual = stampdual))

    # # Parse capacitors # #
    # - Capacitor properties - #
    # name | parent | ID | phases | nominal_voltage | connect_type | capacitorA | capacitorB | capacitorC | pt_phase
    # | delay | upper_limit | lower_limit
    capacitor = []
    capacitor_data = casedata['capacitors']
    if capacitor_data:
        for key, data in capacitor_data.items():
            _id, parent, phases, nom_voltage, connection_type, capacitorA, capacitorB, capacitorC, pt_phase, delay, \
            upper_limit, lower_limit = retrive_attributes(data, 'capacitor', casetype)

            capacitor.append(
                Capacitors(_id, parent, _id, phases, nom_voltage,
                           connection_type, capacitorA, capacitorB, capacitorC,
                           pt_phase, delay, upper_limit, lower_limit))

    # # Parse regulators # #
    # - Regulator properties - #
    # node_index_ | node | name | ID | phases | from_node | to_node | connect_type | reg_type | tapA | tapB | tapC |
    # tapA_change_count | tapB_change_count | tapC_change_count | PT_phase | CT_phase
    regulator = []
    regulator_data = casedata['regulators']
    if regulator_data:
        for key, data in regulator_data.items():
            tap_controls = features[
                "Tap Controls"] if "Tap Controls" in features else None
            tap_controls["Fixed"] if "Tap Controls" in features else True

            _id, phases, from_node, to_node, connection_type, reg_type, tapA, tapB, tapC, pt_phase, \
            ct_phase, connected_transformer, connected_winding, bustype, voltage_setpoint, \
            nom_voltage, band_center, bandwidth, pt_ratio, ct_prim, \
            highstep, lowstep, rated_power, compensator_resistance, compensator_reactance, \
            taps_fixed, aR_min, aR_max, num_tap_steps = retrive_attributes(data, 'regulator', casetype, tap_controls)
            
            # If transformer connected to regulator, create an intermediary connector node.
            if connected_transformer:
                if connected_winding == 1:
                    to_node = to_node + '_RegConnector'
                    flag_connector_exists = next(
                        (obj for obj in node if obj.ID == to_node), None)
                    if not flag_connector_exists:
                        node.append(
                            Nodes(node_index_,
                                  to_node,
                                  to_node,
                                  phases,
                                  nom_voltage,
                                  voltage_setpoint,
                                  connection_type,
                                  False,
                                  bustype=bustype,
                                  voltage_type=voltage_type))
                else:
                    from_node = from_node + '_RegConnector'
                    flag_connector_exists = next(
                        (obj for obj in node if obj.ID == from_node), None)
                    if not flag_connector_exists:
                        node.append(
                            Nodes(node_index_,
                                  from_node,
                                  from_node,
                                  phases,
                                  nom_voltage,
                                  voltage_setpoint,
                                  connection_type,
                                  False,
                                  bustype=bustype,
                                  voltage_type=voltage_type))

            # Get tap control functions if the taps are not fixed
            tap_ctrl_fcn = tap_controls["Function"] if not taps_fixed else None
            tap_ctrl_scheme = tap_controls[
                "Control Scheme"] if not taps_fixed else None
            tap_smoothing_factor = tap_controls[
                "Smoothing Factor"] if not taps_fixed else None

            
            regulator.append(
                Regulators(node_index_, node, _id, _id, phases, from_node,
                          to_node, connection_type, reg_type, nom_voltage,
                          band_center, bandwidth, rated_power, taps_fixed, tapA,
                          tapB, tapC, lowstep, highstep, num_tap_steps, aR_max,
                          aR_min, pt_ratio, pt_phase, ct_prim, ct_phase,
                          connected_transformer, connected_winding,
                          tap_ctrl_scheme, tap_ctrl_fcn, tap_smoothing_factor,
                          compensator_resistance, compensator_reactance, stamp_dual = stampdual))


    # # Switch # #
    # - Switch Properties - #
    # node_index_ | name | ID | from_node | to_node | phases | phase_A_state | phase_B_state | phase_C_state
    switch = []
    switch_data = casedata['switches'] if 'switches' in casedata else None
    if switch_data:
        for key, data in switch_data.items():
            _id, from_node, to_node, phases, phase_A_state, phase_B_state, phase_C_state = retrive_attributes(
                data, 'switch', casetype)
            switch.append(
                Switches(node_index_, _id, _id, from_node, to_node, phases,
                         phase_A_state, phase_B_state, phase_C_state, stamp_dual = stampdual))

    # # Parse fuses # #
    # - Fuse properties - #
    # node_index_ | name | ID | from_node | to_node | phases | phase_A_state | phase_B_state | phase_C_state
    fuse = []
    fuse_data = casedata['fuses'] if 'fuses' in casedata else None
    if fuse_data:
        for key, data in fuse_data.items():
            _id, from_node, to_node, phases, phase_A_state, phase_B_state, phase_C_state = retrive_attributes(
                data, 'fuse')
            fuse.append(
                Fuses(node_index_, _id, _id, from_node, to_node, phases,
                      phase_A_state, phase_B_state, phase_C_state, stamp_dual = stampdual))

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
                Reactors(node, _id, key, from_node, to_node, connection_type,
                         phases, resistance, reactance, stamp_dual = stampdual))

    # DGs nonexistent for now
    dg = []

    return node, slack, load, tplxload, xfmr, ohline, ugline, tplxline, capacitor, regulator, switch, fuse, reactor, dg
