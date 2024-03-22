from django.shortcuts import render, redirect
from django.urls import reverse
from django.views.decorators.csrf import ensure_csrf_cookie
from .models import GLMFile
from .forms import GLMFileForm
import networkx as nx
from django.http import JsonResponse

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
        #uploaded_file = request.FILES['file']  # Make sure 'file' matches the name attribute in your HTML form
        #parsed_data = parse_gml_file(uploaded_file)

        # Redirect to visualization with parsed data
        # This approach uses session to pass data to the next request, adjust as per your requirement
        #request.session['parsed_data'] = parsed_data
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