from __future__ import division, print_function
from django.shortcuts import render, redirect
from django.urls import reverse
from django.views.decorators.csrf import ensure_csrf_cookie
from django.core.files.storage import FileSystemStorage
from .models import GLMFile
from .forms import GLMFileForm
import networkx as nx
from django.http import JsonResponse
from pathlib import Path
import json
import random
import pandas as pd
import pickle
import os
import sys
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
import numpy as np

parent_dir = os.path.dirname(os.getcwd())
sys.path.append(parent_dir + "/classes")
sys.path.append(parent_dir + "/lib")

sys.path.append(parent_dir + "/sugar/classes")
sys.path.append(parent_dir + "/sugar/lib")
sys.path.append(parent_dir + "/forecasting")
sys.path.append(parent_dir + "/forecasting/results")
sys.path.append(parent_dir)

GREEN = "#1cc88a",
RED = "#e74a3b",
GRAY = "#858796",
BLACK = "#2B1B17",
YELLOW = "#f6c23e",
LIGHTBLUE = "#36b9cc",
BLUE = "#4e73df"

import main_forecast

'''for parser'''
from sugar.lib.parser import parser
from sugar.SUGAR3 import main
import lib.read_multi
from default_settings import settings, features
import pyximport
pyximport.install(language_level=3)
# Create settings for parser
infeas_settings_dict = {}
infeas_settings_dict['source type'] = 'current'  # Initial values sent in lib/intialize.py
infeas_settings_dict['obj'] = 'L2'
infeas_settings_dict['obj type'] = 'power'
infeas_settings_dict['obj scalar'] = 1e-6

# VOLTAGE UNBALANCE SETTINGS
voltage_unbalance_settings_dict = {}
# REQUIRED - either True or False
voltage_unbalance_settings_dict['enforce voltage unbalance'] = False

# VOLTAGE BOUND SETTINGS
voltage_bound_settings_dict = {}
# REQUIRED - either True or False
voltage_bound_settings_dict['enforce voltage bounds'] = False

# WARM START SETTINGS
# If performing a warm start to a test, use the dictionary to pass settings
warm_start_settings_dict = {}
# REQUIRED - either True or False
warm_start_settings_dict['initialize'] = False


# feature: battery settings
infeas_settings_dict['battery_node_list'] = [
    {"ID": "B1", "node": "l4", "P_max": 300000, "P_min": 0,
     "Mch": 0.00000005, "Md": 0.00000009, "type": "P", "Bt_prev": 0.1, "C_ch": 1, "C_d": -0.5,
     "single_phase": "A"} # P_max maximum discharge power
]
# node (where to put); ID; Bt_prev: starting charge (use as a percent of maximum charge)
# hardcoded rn (maybe just set the first node to be the battery
SETTINGS = {
    'infeas settings': infeas_settings_dict,
    'Stamp Dual': True,
    'voltage unbalance settings': voltage_unbalance_settings_dict,
    'Warm Start': warm_start_settings_dict,
    'voltage bound settings': voltage_bound_settings_dict,
}
FEATURES = {
    'IBDGs': {
    },
    'Tap Controls': {
        'Fixed': True
    },
    'Current Meas': False,
}



def forecasting_action(request):
    # hours = list(range(1, 49))  # 1 to 48 hours
    # generation = [random.randint(0, 100) for _ in range(1, 49)]
    city = request.session.get('city')
    state = request.session.get('state')
    day_option = request.session.get('day_option')
    casename = request.session.get('casename')
    if day_option == None:
        location = f"{city}, {state}, US"
    else:
        location = None

    # Get Wind Predictions
    wind_pred_df = main_forecast.wind_forecast(location, day_option)
    wind_pred = wind_pred_df['Predictions'].tolist()

    # Get hours once (same for all predictions)
    hours = wind_pred_df['hour']
    hour_labels = hours.apply(lambda x: f"Hour {x}").tolist()

    wind_chart_data = wind_pred

    # Get Solar Predictions
    PV_pred_df = main_forecast.solar_forecast(location, day_option)
    PV_pred = PV_pred_df['Predictions'].tolist()
    solar_chart_data = PV_pred

    # Get Load Predictions
    load_pred_df = main_forecast.load_forecast(location, day_option)
    load_pred = load_pred_df['Predictions'].tolist()
    load_chart_data = load_pred

    context = {
        'wind_chart_labels': json.dumps(hour_labels),
        'wind_chart_data': json.dumps(wind_chart_data),
        'solar_chart_labels': json.dumps(hour_labels),
        'solar_chart_data': json.dumps(solar_chart_data),
        'load_chart_labels': json.dumps(hour_labels),
        'load_chart_data': json.dumps(load_chart_data),
    }
    multicase_name = "temp"
    prices = np.linspace(10,10,24).tolist()
    main_forecast.create_multicase_json(PV_pred, wind_pred, load_pred, prices)
    casename = 'gridlabd/' + casename
    main(casename, multicase_name, settings, features)
     

    return render(request, 'sugarDB/forecasting.html', context)


def visualization_action(request):
    if request.method == 'GET':
        return render(request, 'sugarDB/visualization.html', {'new_upload': "false"})

    return render(request, 'sugarDB/visualization.html', {'new_upload': "false"})


def optimization_action(request):
    base_dir = Path(__file__).resolve().parent.parent.parent
    multi_path = base_dir / 'sugar' / 'multi_outputs.pkl'
    with open(multi_path, 'rb') as file:
        # Load data from the file
        data = pickle.load(file)

    batteryCharge = getRealSumArray(data['P_ch']['B1'][0])
    print("First elements from P_ch:", batteryCharge)

    batteryDischarge = getRealSumArray(data['P_d']['B1'][0])
    print("First elements from P_d:", batteryDischarge)

    batteryStateofCharge = getRealSumArray(data['B']['B1'][0])
    print("First elements from B:", batteryStateofCharge)

    slackGeneration = getRealSumArray(data['P_g']['slack'][0])
    print("First elements from P_g's slack:", slackGeneration)

    renewGeneration = getRealSumArray(data['B_res']['B1'][0])
    print("First elements from P_g's slack:", renewGeneration)

    hours = list(range(len(batteryCharge)))
    batteryCharge_data = pd.DataFrame({
        'hour': hours,
        'batteryCharge': batteryCharge,
    })
    P_ch_chart_labels = batteryCharge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    P_ch_chart_data = batteryCharge_data['batteryCharge'].tolist()

    # battery discharge (P_d)
    hours = list(range(len(batteryDischarge)))
    batteryDischarge_data = pd.DataFrame({
        'hour': hours,
        'batteryDischarge': batteryDischarge,
    })
    P_d_chart_labels = batteryDischarge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    P_d_chart_data = batteryDischarge_data['batteryDischarge'].tolist()

    # battery state of charge (B)
    hours = list(range(len(batteryStateofCharge)))
    batteryStateofCharge_data = pd.DataFrame({
        'hour': hours,
        'batteryStateofCharge': batteryStateofCharge,
    })
    B_chart_labels = batteryStateofCharge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    B_chart_data = batteryStateofCharge_data['batteryStateofCharge'].tolist()

    # slack/diesel generation
    hours = list(range(len(slackGeneration)))
    slackGeneration_data = pd.DataFrame({
        'hour': hours,
        'slackGeneration': slackGeneration,
    })
    sg_chart_labels = slackGeneration_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    sg_chart_data = slackGeneration_data['slackGeneration'].tolist()

    # renewable generation
    hours = list(range(len(renewGeneration)))
    renewGeneration_data = pd.DataFrame({
        'hour': hours,
        'renewGeneration': renewGeneration,
    })
    renew_chart_labels = renewGeneration_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    renew_chart_data = renewGeneration_data['renewGeneration'].tolist()

    context = {
        'P_ch_chart_labels': json.dumps(P_ch_chart_labels),
        'P_ch_chart_data': json.dumps(P_ch_chart_data),
        'P_d_chart_labels': json.dumps(P_d_chart_labels),
        'P_d_chart_data': json.dumps(P_d_chart_data),
        'B_chart_labels': json.dumps(B_chart_labels),
        'B_chart_data': json.dumps(B_chart_data),
        'sg_chart_labels': json.dumps(sg_chart_labels),
        'sg_chart_data': json.dumps(sg_chart_data),
        'renew_chart_labels': json.dumps(renew_chart_labels),
        'renew_chart_data': json.dumps(renew_chart_data),
    }

    return render(request, 'sugarDB/optimization.html', context)

@csrf_exempt
def upload_action(request):
    # Clear the session data when accessing the upload page
    request.session.flush()
    context = {}
    # GET: user navigate here for first time
    if request.method == 'GET':
        context = {'form': GLMFileForm()}
        return render(request, 'sugarDB/upload.html', context)

    # POST: form has been submitted
    form = GLMFileForm(request.POST, request.FILES)
    example_days = request.POST.get('exampleDays')
    if example_days == 'Yes':
        day_option = request.POST.get('dayOptions')
        request.session['city'] = None
        request.session['state'] = None
        request.session['day_option'] = day_option
    else:
        state = request.POST.get('state')
        city = request.POST.get('city')
        request.session['city'] = city
        request.session['state'] = state
        request.session['day_option'] = None

    # print(f"Received: State = {state}, City = {city}")


    # Process the data, save it, send it to another service, etc.
    # print(request.FILES)
    # Validates the form.f
    if form.is_valid():
        # Process the file
        gml_file = request.FILES['file']
        base_dir = Path(__file__).resolve().parent.parent
        testcases_dir = base_dir / 'sugar' / 'testcases' / 'gridlabd'

        filename_base = Path(gml_file.name).stem
        file_dir = testcases_dir / filename_base  # Folder named after the file (without extension)
        #file_dir.mkdir(parents=True, exist_ok=True)  # Create the directory structure
        request.session['casename'] = filename_base

        file_path = file_dir / gml_file.name
        '''with open(file_path, 'wb+') as destination:
            for chunk in gml_file.chunks():
                destination.write(chunk)'''
        try:
            casedata, node_key, node_index_ = parser(str(file_dir), SETTINGS, FEATURES)
            microgrid_data = casedataExtract(casedata)
            #print(casedata.__dict__['node'])
            print(type(casedata))


        except Exception as e:
            print("Parsing failed with exception:", e)
            return redirect('upload')

        microgrid_data_json = json.dumps(microgrid_data)  # Serialize the microgrid data

        return render(request, 'sugarDB/visualization.html', {'microgridData': microgrid_data_json, 'new_upload': "true"})
    else:
        context['form'] = form
        return render(request, 'sugarDB/upload.html', context)


def casedataExtract(casedata):
    microgrid_data = {
        'nodes': [],
        'edges': []
    }
    windloadIDList = []
    PVloadIDList = []
    loadIDList = []
    slackIDList = []

    #---------Load and Nodes------
    for load in casedata.load:
        if load.ID.endswith("_wind"):
            windloadIDList.append(load.ID)
        elif load.ID.endswith("_PV"):
            PVloadIDList.append(load.ID)
        else:
            loadIDList.append(load.ID)

    for slack in casedata.slack:
        slackIDList.append(slack.ID)

    for node in casedata.node:
        if node.ID in loadIDList:
            group = "criticalLoad" # constant power - load._cP * 3 - single phase*3
        elif node.ID in windloadIDList: # nominal power
            group = "windTurbine"
        elif node.ID in PVloadIDList:
            group = "solarPanel"
        elif node.ID in slackIDList:
            group = "generator"
        else:
            group = "Node"

        microgrid_data['nodes'].append({
            'id': node.ID,
            'label': node.name,
            'group': group,
            'value': str(node.Vnom) + "V",
        })

    batteryID = infeas_settings_dict['battery_node_list'][0]['node']
    foundNode = False
    for node in microgrid_data['nodes']:
        if node['id'] == batteryID:
            node['group'] = 'batteryStorage'
            foundNode = True
            break
    if foundNode == False:
        microgrid_data['nodes'].append({
            'id': batteryID,
            'label': infeas_settings_dict['battery_node_list'][0]['node'],
            'group': 'batteryStorage',
            'value': str(infeas_settings_dict['battery_node_list'][0]['P_max']) + "W",
        })

    #----------Edges-------------------
    base_dir = Path(__file__).resolve().parent.parent.parent
    multi_path = base_dir / 'sugar' / 'multi_outputs.pkl'
    with open(multi_path, 'rb') as file:
        # Load data from the file
        data = pickle.load(file)

    for ohline in casedata.ohline:
        MultiOutputs = None
        for line in data['S_line']:
            if line == ohline.name:
                MultiOutputs = getRealSumArray(data['S_line'][line][0])
                # print("Multi output for " + line + " is: ")
                # print(MultiOutputs)

        microgrid_data['edges'].append({
            'from': ohline.from_node,
            'to': ohline.to_node,
            'length': ohline.length,
            'width': 1,
            'label': 'xxx',
            'name': ohline.name,
            'dashes': False,
            'powerFlow': '100',
            'edgetype': 'Overhead lines',
            'multiOutputs': MultiOutputs,
        })

    for ugline in casedata.ugline:
        MultiOutputs = None
        for line in data['S_line']:
            if line == ohline.name:
                MultiOutputs = getRealSumArray(data['S_line'][line][0])
                # print("Multi output for " + line + " is: ")
                # print(MultiOutputs)

        microgrid_data['edges'].append({
            'from': ugline.from_node,
            'to': ugline.to_node,
            'length': ugline.length, #BOUND
            'width': 1,
            'label': 'xxx',
            'name': ugline.name,
            'dashes': True,
            'powerFlow': '100',
            'edgetype': 'Underground lines',
            'multiOutputs': MultiOutputs,
        })


    #-------Special Nodes----------
    for xfmr in casedata.xfmr:
        microgrid_data['edges'].append({
            'from': xfmr.from_node,
            'to': xfmr.to_node,
            'length': 80,
            'width': 1,
            'label': '🌀 transformer 🌀',
            'color': {
                'color': BLUE,
                'highlight': BLUE,
            },
            'font': {
                'color': BLUE,
                'size': 12
            },
            'primaryVoltage': xfmr.primary_voltage,
            'secondaryVoltage': xfmr.secondary_voltage,
            'powerRating': xfmr.power_rating/1000,
            'edgetype': 'Transformer'
        })

    for regulator in casedata.regulator:
        microgrid_data['edges'].append({
            'from': regulator.from_node,
            'to': regulator.to_node,
            'length': 80,
            'width': 1,
            'label': '⚙ regulator ⚙',
            'color': {
                'color': GREEN,
                'highlight': GREEN,
            },
            'font': {
                'color': GREEN,
                'size': 12
            },
            'edgetype': 'Regulator'
        })

    for switch in casedata.switch:
        microgrid_data['edges'].append({
            'from': switch.from_node,
            'to': switch.to_node,
            'length': 80,
            'width': 1,
            'label': '⚡ switch ⚡',
            'color': {
                'color': YELLOW,
                'highlight': YELLOW,
            },
            'font': {
                'color': YELLOW,
                'size': 12
            },
            'color': '#1cc88a',
            'edgetype': 'Switch'
        })



    return microgrid_data



def getRealSumArray(data_line):
    real_parts = np.real(data_line)  # Extract real parts
    abs_real_parts = np.abs(real_parts)  # Compute absolute values
    sum_abs_real_parts = np.sum(abs_real_parts, axis=1).tolist()  # Sum per row
    return sum_abs_real_parts
