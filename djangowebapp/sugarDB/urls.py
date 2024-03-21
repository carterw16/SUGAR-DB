from django.contrib import admin
from django.urls import path, include
from sugarDB import views

urlpatterns = [
    path('', views.visualization_action),
    path('forecasting', views.forecasting_action, name='forecasting'),
    path('upload', views.upload_action, name='upload'),
    path('visualization', views.visualization_action, name='visualization'),
    #path('upload-gml', views.file_upload, name='upload-gml'),
    #path('main', views.main_action, name='main'),
    #path('api/microgrid-data/', views.microgrid_data_view, name='microgrid-data'),
]


