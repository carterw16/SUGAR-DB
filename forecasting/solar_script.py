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
import random_forest, plots

DATA_DIR = os.path.abspath('./data')
results_folder = "results"
# STC solar irradiation for PV modules in W/m^2
GHI_STANDARD = 1000
T_STANDARD = 25
SOLAR_CAPACITY = 400
T_COEFF = -0.003

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
  sunrise = forecast["current"]["sunrise"]
  sunset = forecast["current"]["sunset"]
  sunrise_hour = pd.to_datetime(sunrise, unit='s', utc=True).tz_convert(forecast['timezone']).hour
  sunset_hour = pd.to_datetime(sunset, unit='s', utc=True).tz_convert(forecast['timezone']).hour

  # Initialize lists to store the data
  timestamps = []
  wind_speeds = []
  humidity = []
  temperature = []
  pressure = []
  solar_intensity = []
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
      curr_hour = timestamp.hour

      intensity = np.sin(np.pi * (curr_hour - sunrise_hour) / (sunset_hour - sunrise_hour)) if sunrise_hour <= curr_hour <= sunset_hour else 0
      solar_intensity.append(intensity)
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
  pd.set_option('display.float_format', '{:.4f}'.format)
  return df, pd.DataFrame({"SI": solar_intensity})

def ghi_to_power_factor(GHI, capacity, t_coeff, t_amb):
  GHI_percentage = GHI / GHI_STANDARD
  t_c = GHI * (T_STANDARD/800) + t_amb
  return GHI_percentage - (t_coeff * (t_c - T_STANDARD)) / capacity

def df_ghi_to_power(df):
  p_df = df.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, SOLAR_CAPACITY, T_COEFF, T_STANDARD))
  # pull negatives to 0
  p_df[p_df < 0] = 0
  return p_df

def main():
  filename = 'solarrad_40.44_-79.99_2022.csv'
  X, y = process_training_data(filename)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42, shuffle=True)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  # print(X_test.head())
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")

  # model = lstm_fit(X_train, y_train)
  model = random_forest.train_model(X_train, y_train)
  y_pred = random_forest.predict(model, X_test)
  # Make predictions on the testing data
  # y_pred = lstm_predict(model, X_test)
  power_df = df_ghi_to_power(y_pred)
  # print(power_df.head())
  y_test_power = df_ghi_to_power(y_test)
  y_test_power.columns = ['Power']
  # print(y_test_power.head())
  # power = y_pred.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, 400, -0.03, X_test['Temperature'].mean()))
  # y_test_power = y_test.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, 400, -0.003, X_test['Temperature'].mean()))
  # change columnd header for y_test_power
  # write_predictions(power_df, y_test_power, "solar_predictions.csv")

  # test_metrics = evaluate_model(y_test, y_pred)
  # print("Mean Absolute Error (MAE):", test_metrics['MAE'])
  # print("Mean Squared Error (MSE):", test_metrics['MSE'])
  # print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
  # print("R-squared (R²) Score:", test_metrics['R2'])
  # print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
  # print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
  # write_metrics(test_metrics, 'solar_metrics_rf.txt')

  X_test.reset_index(drop=True, inplace=True)
  # # put hour, X_test, predictions in one dataframe
  plot_df = pd.DataFrame({
    'Hour': X_test['Hour'],
    'Actual': y_test_power['Power'],
    'Predicted': power_df['Predictions'],
    'Temperature': X_test['Temperature']
  })
  # print(plot_df.head())
  # hourly_data = plot_df.groupby('Hour').agg({'Actual':'mean', 'Predicted':'mean'})
  # # keep Hour as a column
  # hourly_data.reset_index(inplace=True)
  # print(hourly_data)

  # Define temperature thresholds for hot and cold days
  hot_threshold = 25  # e.g., days with temperature above 30 degrees Celsius are considered hot
  cold_threshold = 5  # e.g., days with temperature below 10 degrees Celsius are considered cold

  hot_days = plot_df[plot_df['Temperature'] > hot_threshold]
  cold_days = plot_df[plot_df['Temperature'] < cold_threshold]
  hot_hourly = hot_days.groupby('Hour').agg({'Actual':'mean', 'Predicted':'mean'}).reset_index()
  cold_hourly = cold_days.groupby('Hour').agg({'Actual':'mean', 'Predicted':'mean'}).reset_index()
  # Ensure all hours are represented
  all_hours = pd.DataFrame({'Hour': range(24)})
  hot_hourly = pd.merge(all_hours, hot_hourly, on='Hour', how='left').fillna(0)
  cold_hourly = pd.merge(all_hours, cold_hourly, on='Hour', how='left').fillna(0)

  print(hot_hourly)
  print(cold_hourly)
  hourly_data = [hot_hourly, cold_hourly]
  labels = ['Hot Days', 'Cold Days']
  plots.plot_actual_vs_pred(hourly_data, labels)

  # dump(model, "results/solar_model_rf.joblib")
  # # load the model
  # # model = tf.keras.models.load_model("results/solar_model_lstm.keras")
  # model = load("results/solar_model_rf.joblib")
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")
  # forecast_df = format_solar_forecast(forecast)
  # print(forecast_df.head())

  # predictions = random_forest.predict(model, forecast_df)
  # print(predictions.head())

if __name__=="__main__":
    main()
