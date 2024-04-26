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

import os
import sys

parent_dir = os.path.dirname(os.getcwd())
sys.path.append(parent_dir + "/classes")
sys.path.append(parent_dir + "/lib")
sys.path.append(parent_dir + "/forecasting")
sys.path.append(parent_dir + "/forecasting/results")
sys.path.append(parent_dir)


#import main

'''for parser'''
#from sugar.lib.parser import parser
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
     "single_phase": "A"}
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
    '''hours = []
    generations = []
    predictions = main.wind_forecast()
    for predict in predictions['hour']:
        hours.append(predict)
    for predict in predictions['Predictions']:
        generations.append(predict)
    wind_data = pd.DataFrame({
        'hour': hours,
        'generation': generations,
    })
    wind_chart_labels = wind_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    wind_chart_data = wind_data['generation'].tolist()

    predictions = main.solar_forecast()
    for predict in predictions['hour']:
        hours.append(predict)
    for predict in predictions['Predictions']:
        generations.append(predict)
    solar_data = pd.DataFrame({
        'hour': hours,
        'generation': generations,
    })
    solar_chart_labels = solar_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    solar_chart_data = solar_data['generation'].tolist()

    context = {
        'wind_chart_labels': json.dumps(wind_chart_labels),
        'wind_chart_data': json.dumps(wind_chart_data),
        'solar_chart_labels': json.dumps(solar_chart_labels),
        'solar_chart_data': json.dumps(solar_chart_data),
    }'''
    context = {}
    return render(request, 'sugarDB/forecasting.html', context)


def visualization_action(request):
    if request.method == 'GET':
        return render(request, 'sugarDB/visualization.html', {'new_upload': "false"})

    return render(request, 'sugarDB/visualization.html', {'new_upload': "false"})


def optimization_action(request):
    # battery charge (P_ch)
    '''hours = list(range(1, 25))
    batteryCharge = [random.randint(0, 20) for _ in range(1, 25)]
    batteryCharge_data = pd.DataFrame({
        'hour': hours,
        'batteryCharge': batteryCharge,
    })
    P_ch_chart_labels = batteryCharge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    P_ch_chart_data = batteryCharge_data['batteryCharge'].tolist()

    # battery discharge (P_d)
    hours = list(range(1, 25))
    batteryDischarge = [random.randint(0, 20) for _ in range(1, 25)]
    batteryDischarge_data = pd.DataFrame({
        'hour': hours,
        'batteryDischarge': batteryDischarge,
    })
    P_d_chart_labels = batteryDischarge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    P_d_chart_data = batteryDischarge_data['batteryDischarge'].tolist()

    # battery state of charge (B)
    hours = list(range(1, 25))
    batteryStateofCharge = [random.randint(0, 20) for _ in range(1, 25)]
    batteryStateofCharge_data = pd.DataFrame({
        'hour': hours,
        'batteryStateofCharge': batteryStateofCharge,
    })
    B_chart_labels = batteryStateofCharge_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    B_chart_data = batteryStateofCharge_data['batteryStateofCharge'].tolist()

    # slack/diesel generation
    hours = list(range(1, 25))
    slackGeneration = [random.randint(0, 20) for _ in range(1, 25)]
    slackGeneration_data = pd.DataFrame({
        'hour': hours,
        'slackGeneration': slackGeneration,
    })
    sg_chart_labels = slackGeneration_data['hour'].apply(lambda x: f"Hour {x}").tolist()
    sg_chart_data = slackGeneration_data['slackGeneration'].tolist()

    # renewable generation
    renewGeneration = [random.randint(0, 20) for _ in range(1, 25)]
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
'''
    context = {}
    return render(request, 'sugarDB/optimization.html', context)


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
    # print(request.FILES)
    # Validates the form.f
    if form.is_valid():
        # Process the file
        gml_file = request.FILES['file']
        base_dir = Path(__file__).resolve().parent.parent.parent
        testcases_dir = base_dir / 'testcases' / 'gridlabd'

        filename_base = Path(gml_file.name).stem
        file_dir = testcases_dir / filename_base  # Folder named after the file (without extension)
        #file_dir.mkdir(parents=True, exist_ok=True)  # Create the directory structure

        file_path = file_dir / gml_file.name
        '''with open(file_path, 'wb+') as destination:
            for chunk in gml_file.chunks():
                destination.write(chunk)'''
        try:
            #casedata, node_key, node_index_ = parser(str(file_dir), SETTINGS, FEATURES)
            #microgrid_data = casedataExtract(casedata)
            #print(microgrid_data)
            pass

        except Exception as e:
            print("Parsing failed with exception:", e)
            return redirect('upload')

        # microgrid_data_json = json.dumps(microgrid_data)  # Serialize the data
        microgrid_data_json={}
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
            group = "criticalLoad"
        elif node.ID in windloadIDList:
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

    for ohline in casedata.ohline:
        microgrid_data['edges'].append({
            'from': ohline.from_node,
            'to': ohline.to_node,
            'length': ohline.length,
            'width': 4,
            'label': '100 Watts',

            'powerFlow': '100',
            'edgetype': 'Overhead lines'
        })

    for ugline in casedata.ugline:
        microgrid_data['edges'].append({
            'from': ugline.from_node,
            'to': ugline.to_node,
            'length': ugline.length,
            'width': 4,
            'label': '100 Watts',
            'dashes': True,

            'powerFlow': '100',
            'edgetype': 'Underground lines'
        })

    for xfmr in casedata.xfmr:
        microgrid_data['edges'].append({
            'from': xfmr.from_node,
            'to': xfmr.to_node,
            'length': 80,
            'width': 6,
            'label': str() + "V",
            'color': '#4e73df',

            'primaryVoltage': xfmr.primary_voltage,
            'secondaryVoltage': xfmr.secondary_voltage,
            'powerRating': '10',
            'edgetype': 'Transformer'
        })

    for regulator in casedata.regulator:
        microgrid_data['edges'].append({
            'from': regulator.from_node,
            'to': regulator.to_node,
            'length': 80,
            'width': 4,
            'label': 'Regulator',
            'color': '#f6c23e',
            'edgetype': 'Regulator'
        })

    for switch in casedata.switch:
        microgrid_data['edges'].append({
            'from': switch.from_node,
            'to': switch.to_node,
            'length': 80,
            'width': 1,
            'label': 'Switch',
            'color': '#1cc88a',
            'edgetype': 'Switch'
        })



    return microgrid_data


#old
'''def main_action(request):

    # Just display main page form if this is a GET request.
    if request.method == 'GET':
        context = {'form': GLMFileForm()}
        return render(request, 'sugarDB/main.html', context)

    return redirect(reverse('main'))


def file_upload(request):
    context = {}
    # GET: user navigate here for first time
    if request.method == 'GET':
        context = {'form': GLMFileForm()}
        return render(request, 'sugarDB/main.html', context)

    # POST: form has been submitted
    form = GLMFileForm(request.POST, request.FILES)
    # Validates the form.f
    if not form.is_valid():
        context['form'] = form
        return render(request, 'sugarDB/main.html', context)
    else:
        return redirect('main')


# will be called by js ajax request
def microgrid_data_view(request):
    microgrid_data = {
        'nodes': [
            {'id': 1, 'label': 'Controller 1', 'group': 'controller', 'value': 10},
            {'id': 2, 'label': 'Controller 2', 'group': 'controller', 'value': 8},
            {'id': 3, 'label': 'Controller 3', 'group': 'controller', 'value': 6},
            {'id': 4, 'label': 'Controller 3', 'group': 'solarPanel', 'value': 3},
        ],
        'edges': [
            {'from': 1, 'to': 2, 'length': 100, 'width': 6, 'label': '100'},
            {'from': 1, 'to': 3, 'length': 150, 'width': 4, 'label': '150'},
            {'from': 4, 'to': 1, 'length': 150, 'width': 4, 'label': '150'},
            {'from': 4, 'to': 2, 'length': 200, 'width': 4, 'label': '150'},

        ]
    }
    return JsonResponse(microgrid_data)'''