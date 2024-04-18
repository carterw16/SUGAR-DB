import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import write_metrics, write_predictions, evaluate_model, pull_weather_forecast
from lstm import *
import numpy as np
import requests
import urllib.parse
import time
import config
import tensorflow as tf
from joblib import dump, load
import random_forest
# SOLAR_API_KEY = os.getenv("SOLAR_API_KEY")
# EMAIL = "carterw@andrew.cmu.edu"
# BASE_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-5min-download.json?"
# POINTS = [
# '2539138'
# ]
DATA_DIR = os.path.abspath('./data')
results_folder = "results"
# STC solar irradiation for PV modules in W/m^2
GHI_STANDARD = 1000
T_STANDARD = 25
SOLAR_CAPACITY = 400
WEATHER_API_KEY="23e156b4d89df3e0b6c59c6494f7d7cc"
SOLAR_API_KEY = "HCbzILUfWXbR8c6x9YgjSIb6dUWkvGcQL0Yj2Gm2"

def process_training_data(filename):
  data_path = os.path.join(DATA_DIR, filename)
  df = pd.read_csv(data_path)
  df = df.dropna(axis=1, how='any')
  df = df.astype('float64')
  # full_df.set_index('Time', inplace=True)
  # full_df = full_df.apply(pd.to_numeric)
  # full_df = full_df.resample('h').mean().fillna(0)
  # full_df.reset_index(inplace=True)
  # Drop all rows where the minute column is not 0
  df = df[df['Minute'] == 0]
  df['GHI / GHI_standard'] = df['GHI'] / GHI_STANDARD
  # take all columns except GHI and GHI_standard
  features = df.columns.drop(['Year','GHI', 'GHI / GHI_standard'])
  X = df[features]
  # X = df[[Year,Month,Day,Hour,Minute,Wind Speed,Relative Humidity,Temperature,Pressure]]
  y = df["GHI / GHI_standard"]
  return X, y
  # print(df.iloc[-1])
  # cap = 3600
  # full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / cap
  # train_df, test_df = train_test_split(full_df, test_size=0.2)
  # return df

def format_solar_forecast(forecast):
  hourly_weather = forecast["hourly"]

  # Initialize lists to store the data
  timestamps = []
  wind_speeds = []
  humidity = []
  temperature = []
  pressure = []
  # Extract data from hourly_data
  for data_point in hourly_weather:
      # Convert timestamp to datetime object
      timestamp = pd.to_datetime(data_point["dt"], unit="s", utc=True).tz_convert(forecast['timezone'])
      # Append data to lists
      timestamps.append(timestamp)
      wind_speeds.append(data_point["wind_speed"])
      humidity.append(data_point["humidity"])
      temperature.append(data_point["temp"])
      pressure.append(data_point["pressure"])
  # Create DataFrame
  df = pd.DataFrame({
      # "Timestamp": timestamps,
      "Month": [timestamp.month for timestamp in timestamps],
      "Day": [timestamp.day for timestamp in timestamps],
      "Hour": [timestamp.hour for timestamp in timestamps],
      "Minute": [timestamp.minute for timestamp in timestamps],
      "Wind Speed": wind_speeds,
      "Relative Humidity": humidity,
      "Temperature": temperature,
      "Pressure": pressure
  })
  return df

def ghi_to_power_factor(GHI, capacity, t_coeff, t_amb):
  GHI_percentage = GHI / GHI_STANDARD
  t_c = GHI * (T_STANDARD/800) + t_amb
  return GHI_percentage - (t_coeff * (t_c - T_STANDARD)) / capacity

def main():
  filename = 'solarrad_40.44_-79.99_2022.csv'
  X, y = process_training_data(filename)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  print(X_test.head())
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")

  # model = lstm_fit(X_train, y_train)
  model = random_forest.train_model(X_train, y_train)
  y_pred = random_forest.predict(model, X_test)
  # Make predictions on the testing data
  # y_pred = lstm_predict(model, X_test)
  power = y_pred.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, 400, -0.03, X_test['Temperature'].mean()))
  y_test_power = y_test.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, 400, -0.003, X_test['Temperature'].mean()))
  print(y_pred)
  print(y_test)
  # print(power)
  write_predictions(y_pred, y_test, "solar_predictions.csv")

  test_metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", test_metrics['MAE'])
  print("Mean Squared Error (MSE):", test_metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
  print("R-squared (R²) Score:", test_metrics['R2'])
  print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
  write_metrics(test_metrics, 'solar_metrics_rf.txt')

  # model.save(f"results/solar_model_lstm.keras")
  dump(model, "results/solar_model_rf.joblib")
  # load the model
  # model = tf.keras.models.load_model("results/solar_model_lstm.keras")
  model = load("results/solar_model_rf.joblib")
  forecast = pull_weather_forecast("Pittsburgh, PA, US")
  forecast_df = format_solar_forecast(forecast)
  print(forecast_df.head())
  # print(type(forecast_df))

  # predictions = lstm_predict(model, forecast_df)
  predictions = random_forest.predict(model, forecast_df)
  print(predictions.head())

if __name__=="__main__":
    main()
