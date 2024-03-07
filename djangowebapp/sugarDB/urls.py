from django.contrib import admin
from django.urls import path, include
from sugarDB import views

urlpatterns = [
    path('', views.main_action),
    path('upload-gml', views.file_upload, name='upload-gml'),
    path('main', views.main_action, name='main'),
    path('api/microgrid-data/', views.microgrid_data_view, name='microgrid-data'),
]


