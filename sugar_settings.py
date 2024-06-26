####################### METHODS SETTINGS #################################
# MULTIPERIOD SETTINGS
multi_settings_dict = {}
 # REQUIRED - either True or False
multi_settings_dict['run mp'] = True
if multi_settings_dict['run mp']:
    multi_settings_dict['tolerance'] = 1E-3
    multi_settings_dict['verbose'] = True
    multi_settings_dict['max epochs'] = 1


# INFEASIBILITY ANALYSIS SETTINGS
# If performing an infeasibility analysis, use the dictionary to pass analysis settings
infeas_settings_dict = {}
 # REQUIRED - either True or False
infeas_settings_dict['run infeas'] = True
if infeas_settings_dict['run infeas']:
    # OPTIONAL - options are: 'current', 'PQ', 'GB', 'Q', or 'B'; default is 'current'
    infeas_settings_dict['source type'] = 'current' # Initial values sent in lib/intialize.py
    infeas_settings_dict['no print'] = multi_settings_dict['run mp']
    # OPTIONAL - options are: 'L1' or 'L2'; default is 'L2'
    infeas_settings_dict['obj'] = 'L2'
    # OPTIONAL- options are: 'power' or 'current'; default is 'power'
    infeas_settings_dict['obj type'] = 'power'
    # OPTIONAL - scalar value used to scale the objective equation; default is 1; useful when mimizing power or using voltage bounds
    # TODO implement obj scaling value and remove definition in Nodes (self.obj_scalar) and Infeasibility Sources (self.obj_scaling)
    infeas_settings_dict['obj scalar'] = 1e-3
    # OPTIONAL - either True or False; this setting determines if infeasibility sources are added at the Slack bus
    infeas_settings_dict['stamp slack bus'] = True
    # OPTIONAL - default value is 1e-6
    infeas_settings_dict['comp slack tol'] = 1e-1
    # OPTIONAL - either True or False; this setting determines if infeasibility values are printed to the terminal
    infeas_settings_dict['report infeasibility'] = True
    # OPTIONAL - either True or False; this settings determines whether or not to use NetworkX to draw a map of infeasibility in the network
    # TODO: create maps for PQ and GB variables
    infeas_settings_dict['draw infeasibility maps'] = False
    # OPTIONAL - either True or False; this settings is used in mapping - it normalizes infeasibility current
    # TODO: create for PQ and GB variables
    infeas_settings_dict['normalize mapped infeas var'] = False
    # OPTIONAL - either True or False; this setting determines whether infeasibility sources are added at triplex nodes; default value is False
    infeas_settings_dict['triplex sources'] = True # P/Q L1 infeasibility not set up yet
    # OPTIONAL - either True or False; adds an infeasibility source on the neutral
    infeas_settings_dict['neutral infeas source'] = False # Nonzero neutral infeas source not set up yet


# feature: battery settings
infeas_settings_dict['battery_node_list'] = None 
''' [{"ID": "B1", "node":"n3", "P_max":1000, "P_min":0,
                                         "Mch": 1, "Md": 1, "type":"P", "Bt_prev":2000, "C_ch":1, "C_d":-0.5, "single_phase":""}]'''

infeas_settings_dict['slack cost'] = {"w":1}

################################## OPERATIONAL SETTINGS ###########################################
# VOLTAGE BOUND SETTINGS
voltage_bound_settings_dict = {}
# REQUIRED - either True or False
voltage_bound_settings_dict['enforce voltage bounds'] = False
if voltage_bound_settings_dict['enforce voltage bounds']:
    voltage_bound_settings_dict['cs eps tol'] = 1e-3
    voltage_bound_settings_dict['Vmag Bound Max pu'] = 1.05
    voltage_bound_settings_dict['Vmag Bound Min pu'] = 0.95

# VOLTAGE UNBALANCE SETTINGS
voltage_unbalance_settings_dict = {}
# REQUIRED - either True or False
voltage_unbalance_settings_dict['enforce voltage unbalance'] = False
if voltage_unbalance_settings_dict['enforce voltage unbalance']:
    voltage_unbalance_settings_dict['cs eps tol'] = 1e-3
    # The equations are coded s.t. 100 * the upper bound number should equal max percent of voltage unbalance tolerated
    # E.g. .02 as an upper bound corresponds to 2% unbalance
    voltage_unbalance_settings_dict['upper bound'] = 3
    # The lower bound is expected to be zero - if it is not, there will be an error
    voltage_unbalance_settings_dict['lower bound'] = 0

################################## INITIALIZATION SETTINGS ###########################################
# WARM START SETTINGS
# If performing a warm start to a test, use the dictionary to pass settings
warm_start_settings_dict = {}
# REQUIRED - either True or False
warm_start_settings_dict['initialize'] = False
if warm_start_settings_dict['initialize'] == True:
    # REQUIRED: starting from a power flow (0), L1 norm (1), or L2 norm (2) solution?
    warm_start_settings_dict['initial solution type'] = 'L1'
    # REQUIRED: load factor from the previous solution
    warm_start_settings_dict['previous lf'] = 1
    # REQUIRED: objective for new test
    warm_start_settings_dict['new objective'] = 'L1'
    # REQUIRED: file path of the previous solution:
    warm_start_settings_dict['solution file path'] = 'ELWYN_D23802_EXP_VOLTAGES.CSV'
    # REQUIRED: is the file from an external source (e.g. OpenDSS solution)? If so, select True. If it was saved directly from SUGAR, select False
    warm_start_settings_dict["from external simulation"] = True

############################ CALLING SUGAR ###############################################
settings = {
    'Tolerance': 1E-2,
    'Max Iters': 10000,
    'Save File': True,
    'Voltage Limiting':  True,
    'Diode Limiting': True,
    'Load Factor': 1,
    'Homotopy': False,
    'G_homotopy': 400,
    'B_homotopy': 400,
    'Run Type': 'Multi', #'Power Flow', 'Infeas', 'Multi'
    'infeas settings': infeas_settings_dict,
    'multi settings': multi_settings_dict,
    'Stamp Dual': True,
    'Warm Start': warm_start_settings_dict,
    'voltage bound settings': voltage_bound_settings_dict,
    'voltage unbalance settings': voltage_unbalance_settings_dict,
}


features = {
    'IBDGs': {
    },
    'Tap Controls': {
        'Fixed': True
    },
    'Current Meas': True,
}
