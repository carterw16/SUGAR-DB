from django.shortcuts import render, redirect
from django.urls import reverse
from django.views.decorators.csrf import ensure_csrf_cookie
from .models import GLMFile
from .forms import GLMFileForm
import networkx as nx
from django.http import JsonResponse

'''for parser'''
from __future__ import division, print_function
import os
import sys
sys.path.append(os.getcwd() + "/classes")
sys.path.append(os.getcwd() + "/lib")
sys.path.append(os.getcwd())
from SUGAR-DB.lib.parser import parser
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
        uploaded_file = request.FILES['file']  # Make sure 'file' matches the name attribute in your HTML form
        casedata, node_key, node_index_ = parser(uploaded_file, SETTINGS)
        print(casedata.node)
        # Redirect to visualization with parsed data
        # This approach uses session to pass data to the next request, adjust as per your requirement
        request.session['parsed_data'] = casedata.node
        return redirect('visualization')
    else:
        context['form'] = form
        return render(request, 'sugarDB/upload.html', context)


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