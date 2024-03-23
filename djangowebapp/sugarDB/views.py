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
'''for parser'''

import os
import sys
parent_dir = os.path.dirname(os.getcwd())
sys.path.append(parent_dir + "/classes")
sys.path.append(parent_dir + "/lib")
sys.path.append(parent_dir)
from lib.parser import parser
import pyximport
pyximport.install(language_level=3)
# Create settings for parser
infeas_settings_dict = {}
infeas_settings_dict['source type'] = 'current' # Initial values sent in lib/intialize.py
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
                                        {"ID": "B1", "node":"l4", "P_max":300000, "P_min":0,
                                         "Mch": 0.00000005, "Md": 0.00000009, "type":"P", "Bt_prev":0.1, "C_ch":1, "C_d":-0.5, "single_phase":"A"}
]

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
    if request.method == 'GET':
        return render(request, 'sugarDB/forecasting.html')

    return redirect(reverse('forecasting'))


def visualization_action(request):
    if request.method == 'GET':
        return render(request, 'sugarDB/visualization.html')

    return redirect(reverse('visualization'))


def upload_action(request):
    context = {}
    # GET: user navigate here for first time
    if request.method == 'GET':
        context = {'form': GLMFileForm()}
        return render(request, 'sugarDB/upload.html', context)

    # POST: form has been submitted
    form = GLMFileForm(request.POST, request.FILES)
    # Validates the form.f
    if form.is_valid():
        # Process the file
        gml_file = request.FILES['file']  # Make sure 'file' matches the name attribute in your HTML form
        base_dir = Path(__file__).resolve().parent.parent.parent
        testcases_dir = base_dir / 'testcases'

        filename_base = Path(gml_file.name).stem
        file_dir = testcases_dir / filename_base  # Folder named after the file (without extension)
        file_dir.mkdir(parents=True, exist_ok=True)  # Create the directory structure

        file_path = file_dir / gml_file.name
        with open(file_path, 'wb+') as destination:
            for chunk in gml_file.chunks():
                destination.write(chunk)
        microgrid_data = {
            'nodes': [],
            'edges': []
        }

        try:
            print(file_dir)
            # Assuming the parser function requires the full path minus the file extension
            casedata, node_key, node_index_ = parser(parent_dir + '/testcases/gridlabd/13node_ieee_NR_SUGAR', SETTINGS, FEATURES)
            microgrid_data = casedataExtract(casedata)

        except Exception as e:
            print("Parsing failed with exception:", e)
            return redirect('upload')

        microgrid_data_json = json.dumps(microgrid_data)  # Serialize the data
        return render(request, 'sugarDB/visualization.html', {'microgridData': microgrid_data_json})
    else:
        context['form'] = form
        return render(request, 'sugarDB/upload.html', context)


def casedataExtract(casedata):
    microgrid_data = {
        'nodes': [],
        'edges': []
    }

    for node in casedata.node:
        group = "controller"
        if node.name[0] == "L":
            group = "criticalLoad"

        microgrid_data['nodes'].append({
            'id': node.ID,
            'label': node.name,
            'group': group,
            'value': node.Vnom
        })


    for ohline in casedata.ohline:
        microgrid_data['edges'].append({
            'from': ohline.from_node,
            'to': ohline.to_node,
            'length': ohline.length,
            'width': 4,
            'label': ohline.freq
        })

    for ugline in casedata.ugline:
        microgrid_data['edges'].append({
            'from': ugline.from_node,
            'to': ugline.to_node,
            'length': ugline.length,
            'width': 1,
            'label': ugline.freq
        })

    for xfmr in casedata.xfmr:
        microgrid_data['edges'].append({
            'from': xfmr.from_node,
            'to': xfmr.to_node,
            'length': 200,
            'width': 6,
            'label': xfmr.primary_voltage
        })

    for regulator in casedata.regulator:
        microgrid_data['edges'].append({
            'from': regulator.from_node,
            'to': regulator.to_node,
            'length': 100,
            'width': 4,
            'label': regulator.Vmax
        })

    # fuse?tplxline?

    for switch in casedata.switch:
        microgrid_data['edges'].append({
            'from': switch.from_node,
            'to': switch.to_node,
            'length': 300,
            'width': 1,
            'label': switch.phases
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