import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import evaluate_model
import numpy as np
import requests
from datetime import datetime
import urllib.parse
import time
from wind_script import evaluate_model, write_metrics, write_predictions, pull_weather_forecast
import config
import tensorflow as tf
from lstm import *

WEATHER_API_KEY = config.WEATHER_API_KEY
DATA_DIR = os.path.abspath('./data')
MEAN_LOAD = 421

def process_training_data(load_file, weather_file):
  load_path = os.path.join(DATA_DIR, load_file)
  weather_path = os.path.join(DATA_DIR, weather_file)
  load_df = pd.read_csv(load_path)

  load_df['Time'] = pd.to_datetime(load_df['Time'])
  load_df = load_df[['Time','Aggregate']]
  load_df = load_df.dropna(axis=0, how='any')
  load_df.set_index('Time', inplace=True)
  load_df = load_df.resample('h').mean().fillna(0)
  load_df.reset_index(inplace=True)

  # print average of aggregate column
  # print("Average load: ", load_df['Aggregate'].mean())
  # print(df['Aggregate'])

  weather_df = pd.read_csv(weather_path)
  weather_df = weather_df[['dt', 'temp', 'pressure', 'humidity', 'clouds_all']]
  weather_df['dt'] = pd.to_datetime(weather_df['dt'], unit='s')
  weather_df=weather_df.rename(columns={'dt': 'Time'})

  df = load_df.merge(weather_df, on='Time')
  df = load_df
  # max_load = df['Aggregate'].max()
  # print(max_load)
  df['Aggregate'] = df['Aggregate'] / MEAN_LOAD
  # df['month'] = df['Time'].dt.month
  # df['day'] = df['Time'].dt.dayofweek
  df['hour'] = df['Time'].dt.hour
  X = df.drop(columns=['Aggregate', 'Time'])
  y = df['Aggregate']
  return X,y

def format_load_forecast(forecast):
  hourly_weather = forecast["hourly"]
  timestamps = []
  # Extract data from hourly_data
  for data_point in hourly_weather:
      timestamp = pd.to_datetime(data_point["dt"], unit="s", utc=True).tz_convert(forecast['timezone'])
      timestamps.append(timestamp)
  df = pd.DataFrame({
      # "Timestamp": timestamps,
      # "month": [timestamp.month for timestamp in timestamps],
      # "dayofweek": [timestamp.dayofweek for timestamp in timestamps],
      "hour": [timestamp.hour for timestamp in timestamps],
  })
  return df

def load_to_power(df):
  return df * MEAN_LOAD

def main():
  load_file = 'UK_house_loads/House_1.csv'
  weather_file = 'openweather_hourly_2013_2023/loughborough.csv'
  X,y = process_training_data(load_file, weather_file)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")
  print(X_train.head())

  # # linear_regression(X_train, X_test, y_train, y_test)
  model = lstm_fit(X_train, y_train)
  # Make predictions on the testing data
  y_pred = lstm_predict(model, X_test)
  write_predictions(y_pred, y_test, "load_predictions.csv")

  test_metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", test_metrics['MAE'])
  print("Mean Squared Error (MSE):", test_metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
  print("R-squared (R²) Score:", test_metrics['R2'])
  print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
  write_metrics(test_metrics, 'load_metrics_lstm.txt')

  model.save(f"results/load_model_lstm.keras")

  model = tf.keras.models.load_model("results/load_model_lstm.keras")
  forecast = pull_weather_forecast("Pittsburgh, PA, US")
  load_df = format_load_forecast(forecast)
  load_predictions = lstm_predict(model, load_df)
  print(load_predictions)
  load_predictions = load_to_power(load_predictions)
  print(load_predictions)

if __name__=="__main__":
  main()