import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import evaluate_model
from sklearn.preprocessing import StandardScaler, RobustScaler
import numpy as np
import requests
from datetime import datetime
import urllib.parse
import time
from wind_script import evaluate_model, write_metrics, write_predictions, pull_weather_forecast, linear_regression_fit, linear_regression_predict
import config
import tensorflow as tf
from lstm import *
import plots, gradient_boosting, random_forest
from joblib import dump, load

DATA_DIR = os.path.abspath('./data')
# MEAN_LOAD = 421
# MAX_LOAD = 2895
MAX_LOAD = 1000

def preprocess_data(load_df, weather_df, scaler_x, scaler_y):
    load_df['Time'] = pd.to_datetime(load_df['Time'])
    load_df.set_index('Time', inplace=True)
    df = load_df
    weather_df['Time'] = pd.to_datetime(weather_df['dt'], unit='s')
    weather_df.set_index('Time', inplace=True)
    # Drop all columns from weather_df except the index, feels_like and humidity
    weather_df = weather_df[['feels_like', 'humidity']]
    # Merge datasets
    df = pd.merge(load_df, weather_df, left_index=True, right_index=True, how='inner')
    # Drop any rows with NaN values that resulted from resampling
    df.dropna(inplace=True)

    print(df.head())

    # Extract temporal features from the index
    df['hour'] = df.index.hour
    df['day_of_week'] = df.index.dayofweek
    df['month'] = df.index.month
    # df['day_of_month'] = df.index.day

    df['AverageLoad'] = df['AverageLoad'] / scaler_y
    df.reset_index(inplace=True)
    df.drop(columns=['Time'], inplace=True)

    # Separate features and target
    X = df.drop('AverageLoad', axis=1)
    y = df['AverageLoad']

    return X, y


def format_load_forecast(forecast):
  hourly_weather = forecast["hourly"]
  timestamps = []
  feels_like = []
  humidity = []
  # Extract data from hourly_data
  for data_point in hourly_weather:
      timestamp = pd.to_datetime(data_point["dt"], unit="s", utc=True).tz_convert(forecast['timezone'])
      timestamps.append(timestamp)
      feels_like.append(data_point["feels_like"])
      humidity.append(data_point["humidity"])
  df = pd.DataFrame({
      # "Time": timestamps,
      "hour": [timestamp.hour for timestamp in timestamps],
      "day_of_week": [timestamp.dayofweek for timestamp in timestamps],
      "month": [timestamp.month for timestamp in timestamps],
      "feels_like": feels_like,
      "humidity": humidity
  })
  # df.set_index('Time', inplace=True)
  # features = ['hour']
  # print(df.head())
  # df[features] = scaler_x.fit_transform(df[features].values)
  return df

def collect_load_data():
  # get a list of paths to the csv files in data/UK_house_loads
  houses = DATA_DIR + '/UK_house_loads'
  file_paths = [os.path.join(houses, file) for file in os.listdir(houses) if file.endswith('.csv')]
  print(file_paths)

  # List to store the hourly average DataFrames
  hourly_averages = []

  for file_path in file_paths:
      # Read the CSV file
      df = pd.read_csv(file_path, parse_dates=['Time'], index_col='Time')

      # Resample to hourly and calculate the mean of the "Aggregate" column
      hourly_avg = df['Aggregate'].resample('H').mean()

      # Store the hourly averages DataFrame
      hourly_averages.append(hourly_avg)

  # Combine the hourly averages into a single DataFrame
  combined_hourly_avg = pd.concat(hourly_averages, axis=1)

  # If you want the average across all houses for each hour
  average_across_houses = combined_hourly_avg.mean(axis=1)

  # Save the combined hourly averages to a new CSV file
  average_across_houses.to_csv(DATA_DIR + '/average_loads.csv', header=['AverageLoad'])

def plot_load_results(X_actual, y_actual, predictions):
  # put hour, X_test, predictions in one dataframe
  plot_df = pd.DataFrame({
    'Hour': X_actual['hour'],
    'Actual': y_actual['AverageLoad'],
    'Predicted': predictions['Predictions'],
    # 'Temperature': load_df['Temperature']
  })
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

def main():
  load_file = 'average_loads.csv'
  weather_file = 'openweather_hourly_2013_2023/loughborough.csv'
  load_path = os.path.join(DATA_DIR, load_file)
  weather_path = os.path.join(DATA_DIR, weather_file)
  load_df = pd.read_csv(load_path)
  weather_df = pd.read_csv(weather_path)
  scaler_x = MinMaxScaler()
  X, y = preprocess_data(load_df, weather_df, scaler_x, MAX_LOAD)
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, shuffle=True)
  y_test = pd.DataFrame(y_test)
  # print(X_test.head())
  y_test.reset_index(drop=True, inplace=True)

  model = random_forest.train_model(X_train, y_train)

  predictions = random_forest.predict(model, X_test)

  test_metrics = evaluate_model(y_test, predictions)

  write_predictions(predictions, y_test, "load_predictions.csv")
  write_metrics(test_metrics, 'load_metrics_rf.txt')
  dump(model, "results/load_model_rf.joblib")

  X_test.reset_index(drop=True, inplace=True)
  # # put hour, X_test, predictions in one dataframe
  plot_df = pd.DataFrame({
    'Hour': X_test['hour'],
    'Actual': y_test['AverageLoad'],
    'Predicted': predictions['Predictions'],
    'Temperature': X_test['feels_like']
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

  model = load("results/load_model_rf.joblib")

  forecast = pull_weather_forecast("Pittsburgh, PA, US")
  load_df = format_load_forecast(forecast)
  load_predictions = random_forest.predict(model, load_df)
  load_predictions['hour'] = load_df.hour
  print(load_predictions)

if __name__=="__main__":
  main()
