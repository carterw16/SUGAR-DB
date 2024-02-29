from django.shortcuts import render
from django.views.decorators.csrf import ensure_csrf_cookie


def main_action(request):

    # Just display main page form if this is a GET request.
    if request.method == 'GET':
        return render(request, 'sugarDB/main.html')

    return render(request, 'surgarDB/main.html')
