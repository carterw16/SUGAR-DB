"""Saves the results from a simulation.

  Author(s): Naeem Turner-Bandele
  Created Date: 10-02-2019
  Updated Date: 10-18-2020
  Email: nturnerb@cmu.edu
  Status: Development

  Saves simulation results into CSVs. Currently, this function is designed to show voltages, currents,
  regulator-specific outputs, and IBDG-specific outputs.

"""

import copy
import csv
import os
from pathlib import Path

import numpy as np
from natsort import natsorted

from classes.Nodes import Nodes

def save_output(output_dir, case, node, outputs, iterationCount, 
                simulation_stats, settings, enable_IBDGs, lf, stamp_dual):
    casename = None
    if 'gridlabd' in case:
        casename = case.replace('gridlabd/', '')
    elif 'opendss' in case:
        casename = case.replace('opendss/', '')
    elif 'cyme' in case:
        casename = case.replace('cyme/', '')

    output_subdir = output_dir + os.path.sep + casename + os.path.sep
    path_to_output = output_subdir + os.path.basename(casename)
    Path(output_subdir).mkdir(parents=True, exist_ok=True)

    load = outputs[0]
    ibdg = outputs[1]
    regulator = outputs[2]
    current_meas = outputs[3]
    transformer = outputs[4]
    fuse = outputs[5]
    switch = outputs[6]

    summary_path = '_' + str(lf) + '_summary_SUGAR.csv'
    if enable_IBDGs:
        ibdg_type = str(type(outputs[1][0])).replace("<class 'classes.IBDGs.",
                                                     "")
        ibdg_type = ibdg_type.replace("'>", "")
        volts_path = '_' + ibdg_type + '_' + str(lf) + '_voltages_SUGAR.csv'
        currents_path = '_' + ibdg_type + '_' + str(lf) + '_currents_SUGAR.csv'
        '_' + ibdg_type + '_' + str(lf) + '_loads_SUGAR.csv'
    else:
        volts_path = '_' + str(lf) + '_voltages_SUGAR.csv'
        currents_path = '_' + str(lf) + '_currents_SUGAR.csv'
        '_' + str(lf) + '_loads_SUGAR.csv'
    
    if stamp_dual:
        dual_analysis_ABC = outputs[8]
        dual_analysis_tplx = outputs[9]
        infeas_path = '_' + str(lf) + '_infeas_SUGAR.csv'
        infeas_path_tplx = '_' + str(lf) + '_infeas_SUGAR_tplx.csv'
    

    # Simulation Stats
    reader_time = np.around(simulation_stats[0], 8)
    modeler_time = np.around(simulation_stats[1], 8)
    parser_time = np.around(simulation_stats[2], 8)
    sim_total_time = np.around(simulation_stats[3], 5)

    # Total Load
    P_total = 0
    Q_total = 0
    for ele in load:
        P = np.sum(np.asarray(ele._cP))
        Q = np.sum(np.asarray(ele._cQ))
        P_total += P
        Q_total += Q
    total_load = np.abs(complex(P_total, Q_total))

    # Total Solar Gen
    P_total = 0
    Q_total = 0
    for ele in ibdg:
        P = float(ele.P3)
        Q = float(ele.Q3_final)
        P_total += P
        Q_total += Q
    total_solar = np.abs(complex(P_total, Q_total))

    # Total Nodes
    total_nodes = Nodes.voltage_index.__len__() + Nodes.Q_index.__len__() + Nodes.Tap_index.__len__()

    with open(path_to_output + summary_path, 'w', newline='\n') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(['SIMULATION PROFILE AND SUMMARY'])
        csvwriter.writerow(['Iters', iterationCount])
        csvwriter.writerow(['Reader Time (s)', reader_time])
        csvwriter.writerow(['Modeler Time (s)', modeler_time])
        csvwriter.writerow(['Parse Time (s)', parser_time])
        csvwriter.writerow(['Solve Time (s)', sim_total_time])
        csvwriter.writerow(['Homotopy Enabled', settings['Homotopy']])
        csvwriter.writerow(['Voltage Limiting', settings['Voltage Limiting']])
        csvwriter.writerow(['Number of Buses', node.__len__()])
        csvwriter.writerow(['Number of Nodes', total_nodes])
        csvwriter.writerow(['Total Generation', total_solar])
        csvwriter.writerow(['Total Load', total_load])

    with open(path_to_output + volts_path, 'w', newline='\n') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        csvwriter.writerow(['LINE-GROUND VOLTAGES BY BUS & PHASE'])
        _node_temp = copy.copy(node)
        csvwriter.writerow([
            'Bus', 'BasekV', 'phase1', 'Magnitude1', 'Angle1', 'pu1',
            'phase2', 'Magnitude2', 'Angle2', 'pu2', 'phase3',
            'Magnitude3', 'Angle3', 'pu3', ' Unbalance (%)'
        ])

        for ele in _node_temp:
            _temp_var = ele.name.split('_')
            # For Taxonomical cases
            if len(_temp_var) == 3:
                if len(_temp_var[2]) == 1:
                    _temp_var[2] = '0' + '0' + '0' + _temp_var[2]
                elif len(_temp_var[2]) == 2:
                    _temp_var[2] = '0' + '0' + _temp_var[2]
                elif len(_temp_var[2]) == 3:
                    _temp_var[2] = '0' + _temp_var[2]
                _temp = "".join(_temp_var)
                ele.name_temp = _temp
            else:
                ele.name_temp = ele.name

        SortedNames = natsorted(enumerate(node),
                                key=lambda y: y[1].name_temp.lower())

        n_names = len(SortedNames)
        for ele in range(n_names):
            node_ele = node[SortedNames[ele][0]]
            if not node_ele.isTriplex:
                if node_ele.ID.find('_RegConnector') == -1:
                    csvwriter.writerow([
                        node_ele.name, node_ele.Vbase_LL, 'A', node_ele.Vln[0],
                        node_ele.V_ang[0], node_ele.Vln_pu[0], 'B',
                        node_ele.Vln[1], node_ele.V_ang[1], node_ele.Vln_pu[1],
                        'C', node_ele.Vln[2], node_ele.V_ang[2],
                        node_ele.Vln_pu[2], node_ele.V_unb
                    ])
            else:
                csvwriter.writerow([
                    node_ele.name, node_ele.Vbase_LL, 1,  node_ele.Vln[0],
                    node_ele.V_ang[0], node_ele.Vln_pu[0], 2, node_ele.Vln[1],
                    node_ele.V_ang[1], node_ele.Vln_pu[1], 0, 0.0, 0.0, 0.0,
                    node_ele.V_unb
                ])

    if ibdg:
        ibdg_path = '_' + str(lf) + '_ibdgs_SUGAR.csv'
        with open(path_to_output + ibdg_path, 'w', newline='\n') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow([
                'der_name', 'sensor', 'type', 'currA_mag', 'currA_ang',
                'currB_mag', 'currB_ang', 'currC_mag', 'currC_ang', 'Qmin',
                'Q3', 'Qmax', '%Q'
            ])
            for ele in ibdg:
                csvwriter.writerow([
                    ele.name, ele.Sensor, ibdg_type, ele.Ia_mag, ele.Ia_ang,
                    ele.Ib_mag, ele.Ib_ang, ele.Ic_mag, ele.Ic_ang, ele.Qmin,
                    ele.Q3_final, ele.Qmax, ele.perc_Q
                ])

    if regulator:
        regulator_path = '_' + str(lf) + '_regulators_SUGAR.csv'
        with open(path_to_output + regulator_path, 'w',
                  newline='\n') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(['REGULATORS'])
            csvwriter.writerow([
                'Name', 'Tap A', 'Tap B', 'Tap C', 'aR A', 'aR B', 'aR C',
                'Ia Magnitude', 'Ia Angle',
                'Ib Magnitude', 'Ib Angle',
                'Ic Magnitude', 'Ic Angle'
            ])
            for ele in regulator:
                csvwriter.writerow([
                    ele.name, ele.tapA, ele.tapB, ele.tapC, ele.aR_A, ele.aR_B,
                    ele.aR_C, ele.Ia_mag, ele.Ia_ang, ele.Ib_mag, ele.Ib_ang,
                    ele.Ic_mag, ele.Ic_ang
                ])

    # Save Measured Line Currents
    with open(path_to_output + currents_path, 'w', newline='\n') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        if current_meas:
            csvwriter.writerow(['LINE CURRENTS BY CONNECTION & PHASE'])
            csvwriter.writerow([
                'Connection', 'From', 'To', 'Magnitude1', 'Angle1',
                'Magnitude2', 'Angle2', 'Magnitude3', 'Angle3'
            ])
            for ele in current_meas:
                if not ele.isTriplex:
                    csvwriter.writerow([
                        ele.name, ele.from_node, ele.to_node, ele.Ia_mag,
                        ele.Ia_ang, ele.Ib_mag, ele.Ib_ang, ele.Ic_mag,
                        ele.Ic_ang
                    ])
                else:
                    csvwriter.writerow([
                        ele.name, ele.from_node, ele.to_node, ele.I1_mag,
                        ele.I1_ang, ele.I2_mag, ele.I2_ang, 0.0, 0.0
                    ])

            for ele in transformer:
                if ele.connect_type != 5:
                    csvwriter.writerow([
                        ele.name, ele.from_node, ele.to_node, ele.Ia_mag,
                        ele.Ia_ang, ele.Ib_mag, ele.Ib_ang, ele.Ic_mag,
                        ele.Ic_ang
                    ])
                else:
                    csvwriter.writerow([
                        ele.name, ele.from_node, ele.to_node, ele.I1_mag,
                        ele.I1_ang, ele.I2_mag, ele.I2_ang, 0.0, 0.0
                    ])

            for ele in fuse:
                csvwriter.writerow([
                    ele.name, ele.from_node, ele.to_node, ele.Ia_mag,
                    ele.Ia_ang, ele.Ib_mag, ele.Ib_ang, ele.Ic_mag, ele.Ic_ang
                ])

            for ele in switch:
                csvwriter.writerow([
                    ele.name, ele.from_node, ele.to_node, ele.Ia_mag,
                    ele.Ia_ang, ele.Ib_mag, ele.Ib_ang, ele.Ic_mag, ele.Ic_ang
                ])
    
    if stamp_dual:
        with open(path_to_output + infeas_path, 'w', newline='\n') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(['Infeasibility Currents'])
            csvwriter.writerow([dual_analysis_ABC.columns])
            if 'Infeas P (A)' in dual_analysis_ABC.columns and 'Infeas Q (A)' in dual_analysis_ABC.columns:
                for i in range(1,len(dual_analysis_ABC)):
                    csvwriter.writerow([
                        dual_analysis_ABC.index[i], 
                        dual_analysis_ABC['Node Name'][i],
                        dual_analysis_ABC['Infeas Ir (A)'][i],
                        dual_analysis_ABC['Infeas Ii (A)'][i],
                        dual_analysis_ABC['Infeas Ir (B)'][i],
                        dual_analysis_ABC['Infeas Ii (B)'][i],
                        dual_analysis_ABC['Infeas Ir (C)'][i],
                        dual_analysis_ABC['Infeas Ii (C)'][i],
                        dual_analysis_ABC['Infeas P (A)'][i],
                        dual_analysis_ABC['Infeas Q (A)'][i],
                        dual_analysis_ABC['Infeas P (B)'][i],
                        dual_analysis_ABC['Infeas Q (B)'][i],
                        dual_analysis_ABC['Infeas P (C)'][i],
                        dual_analysis_ABC['Infeas Q (C)'][i],                       
                    ])   
            elif 'Infeas P (A)' not in dual_analysis_ABC.columns and 'Infeas Q (A)' in dual_analysis_ABC.columns:
                for i in range(1,len(dual_analysis_ABC)):    
                    csvwriter.writerow([
                        dual_analysis_ABC.index[i], 
                        dual_analysis_ABC['Node Name'][i],
                        dual_analysis_ABC['Infeas Ir (A)'][i],
                        dual_analysis_ABC['Infeas Ii (A)'][i],
                        dual_analysis_ABC['Infeas Ir (B)'][i],
                        dual_analysis_ABC['Infeas Ii (B)'][i],
                        dual_analysis_ABC['Infeas Ir (C)'][i],
                        dual_analysis_ABC['Infeas Ii (C)'][i],
                        dual_analysis_ABC['Infeas Q (A)'][i],
                        dual_analysis_ABC['Infeas Q (B)'][i],
                        dual_analysis_ABC['Infeas Q (C)'][i],                       
                    ])  
            elif 'Infeas G (A)' in dual_analysis_ABC.columns and 'Infeas B (A)' in dual_analysis_ABC.columns:
                for i in range(1,len(dual_analysis_ABC)):    
                    csvwriter.writerow([
                        dual_analysis_ABC.index[i], 
                        dual_analysis_ABC['Node Name'][i],
                        dual_analysis_ABC['Infeas Ir (A)'][i],
                        dual_analysis_ABC['Infeas Ii (A)'][i],
                        dual_analysis_ABC['Infeas Ir (B)'][i],
                        dual_analysis_ABC['Infeas Ii (B)'][i],
                        dual_analysis_ABC['Infeas Ir (C)'][i],
                        dual_analysis_ABC['Infeas Ii (C)'][i],
                        dual_analysis_ABC['Infeas G (A)'][i],
                        dual_analysis_ABC['Infeas B (A)'][i],
                        dual_analysis_ABC['Infeas G (B)'][i],
                        dual_analysis_ABC['Infeas B (B)'][i],
                        dual_analysis_ABC['Infeas G (C)'][i],
                        dual_analysis_ABC['Infeas B (C)'][i],                       
                    ])              
            elif 'Infeas G (A)' not in dual_analysis_ABC.columns and 'Infeas B (A)' in dual_analysis_ABC.columns:
                for i in range(1,len(dual_analysis_ABC)):
                    csvwriter.writerow([
                        dual_analysis_ABC.index[i], 
                        dual_analysis_ABC['Node Name'][i],
                        dual_analysis_ABC['Infeas Ir (A)'][i],
                        dual_analysis_ABC['Infeas Ii (A)'][i],
                        dual_analysis_ABC['Infeas Ir (B)'][i],
                        dual_analysis_ABC['Infeas Ii (B)'][i],
                        dual_analysis_ABC['Infeas Ir (C)'][i],
                        dual_analysis_ABC['Infeas Ii (C)'][i],
                        dual_analysis_ABC['Infeas B (A)'][i],
                        dual_analysis_ABC['Infeas B (B)'][i],
                        dual_analysis_ABC['Infeas B (C)'][i],                       
                    ])
            else:
                for i in range(1,len(dual_analysis_ABC)):
                    csvwriter.writerow([
                        dual_analysis_ABC.index[i], 
                        dual_analysis_ABC['Node Name'][i],
                        dual_analysis_ABC['Infeas Ir (A)'][i],
                        dual_analysis_ABC['Infeas Ii (A)'][i],
                        dual_analysis_ABC['Infeas Ir (B)'][i],
                        dual_analysis_ABC['Infeas Ii (B)'][i],
                        dual_analysis_ABC['Infeas Ir (C)'][i],
                        dual_analysis_ABC['Infeas Ii (C)'][i]
                    ])   
        if dual_analysis_tplx is not None:
            with open(path_to_output + infeas_path_tplx, 'w', newline='\n') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',')
                csvwriter.writerow(['Infeasibility Currents'])
                csvwriter.writerow([
                    'Node Name', 'Ifr, 1', 'Ifi, 1', 'Ifr,2', 'Ifi,2', 'Ifr,N', 'Ifi,N'
                ])
                for i in range(1,len(dual_analysis_tplx)):
                    csvwriter.writerow([
                        dual_analysis_tplx.index[i], 
                        dual_analysis_tplx['Node Name'][i],
                        dual_analysis_tplx['if1,r Value'][i],
                        dual_analysis_tplx['if1,i Value'][i],
                        dual_analysis_tplx['if2,r Value'][i],
                        dual_analysis_tplx['if2,i Value'][i],
                        dual_analysis_tplx['ifN,r Value'][i],
                        dual_analysis_tplx['ifN,i Value'][i]
                    ])               

