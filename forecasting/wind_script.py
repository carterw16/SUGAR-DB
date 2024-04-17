import pandas as pd
import os, datetime, requests
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from keras.models import Sequential
from keras.layers import LSTM, Dense, Input
import tensorflow as tf
import config
from config import *
from lstm import *
# from automl import *
import random_forest
from joblib import dump, load

WEATHER_API_KEY = "23e156b4d89df3e0b6c59c6494f7d7cc"
DATA_DIR = os.path.abspath('./data')
results_folder = "results"
WIND_CAPACITY = 3600

def process_training_data(filename):
  data_path = os.path.join(DATA_DIR, filename)
  full_df = pd.read_csv(data_path)
  # cap = 3600
  full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / WIND_CAPACITY
  full_df['Time'] = pd.to_datetime(full_df['Date/Time'], format='%d %m %Y %H:%M')
  # drop the Date/time column from full_df
  full_df.drop(columns=['Date/Time'], inplace=True)
  full_df.set_index('Time', inplace=True)
  full_df = full_df.apply(pd.to_numeric)
  full_df = full_df.resample('h').mean().fillna(0)
  full_df.reset_index(inplace=True)
  # clip wind speeds at 2 decimal places
  full_df['Wind Speed (m/s)'] = full_df['Wind Speed (m/s)'].apply(lambda x: round(x, 2))
  full_df['month'] = full_df['Time'].dt.month
  full_df['day'] = full_df['Time'].dt.dayofweek
  full_df['hour'] = full_df['Time'].dt.hour
  full_df = full_df[['month', 'day', 'hour', 'Wind Speed (m/s)', 'LV ActivePower (%)']]
  X = full_df[["month","day","hour","Wind Speed (m/s)"]]
  y = full_df["LV ActivePower (%)"]
  # train_df, test_df = train_test_split(full_df, test_size=0.2, shuffle=False)
  # return train_df, test_df
  return X, y

def format_wind_forecast(forecast):
  hourly_weather = forecast["hourly"]

  # Initialize lists to store the data
  timestamps = []
  wind_speeds = []
  wind_degrees = []

  # Extract data from hourly_data
  for data_point in hourly_weather:
      # Convert timestamp to datetime object
      timestamp = pd.to_datetime(data_point["dt"], unit="s", utc=True).tz_convert(forecast['timezone'])
      # Append data to lists
      timestamps.append(timestamp)
      wind_speeds.append(data_point["wind_speed"])
      wind_degrees.append(data_point["wind_deg"])

  # Create DataFrame
  df = pd.DataFrame({
      # "Timestamp": timestamps,
      "month": [timestamp.month for timestamp in timestamps],
      "day": [timestamp.dayofweek for timestamp in timestamps],
      "hour": [timestamp.hour for timestamp in timestamps],
      "Wind Speed (m/s)": wind_speeds,
      # "Wind Direction (°)": wind_degrees
  })
  return df

def pull_weather_forecast(city):
  cnt = 24  # Number of hourly data points to retrieve per day

  geocoding_url = "http://api.openweathermap.org/geo/1.0/direct"

  # Parameters for the Geocoding API request
  params = {
      "q": city,
      "limit": 1,
      "appid": WEATHER_API_KEY
  }

  # Send the request
  response = requests.get(geocoding_url, params=params)

  # Check if request was successful
  if response.status_code == 200:
      # Parse the JSON response
      data = response.json()
      # Extract latitude and longitude
      lat = data[0]["lat"]
      lon = data[0]["lon"]
  else:
      print("Failed to fetch geocoding data:", response.status_code)
      return

  onecall_url = "https://api.openweathermap.org/data/3.0/onecall"

  # Parameters for the One Call API request
  params = {
      "lat": lat,
      "lon": lon,
      "exclude": "minutely,daily",
      "units": "metric",
      "appid": WEATHER_API_KEY
  }

  # Send the request
  response = requests.get(onecall_url, params=params)

  # Check if request was successful
  if response.status_code == 200:
    data = response.json()
    return data
  else:
    print("Failed to fetch weather data:", response.status_code)
    return

def evaluate_model(ground_truth, predictions):
    range = np.ptp(ground_truth)
    metrics = {}
    metrics['MAE'] = mean_absolute_error(ground_truth, predictions)
    metrics['MSE'] = mean_squared_error(ground_truth, predictions)
    metrics['RMSE'] = np.sqrt(metrics['MSE'])
    metrics['R2'] = r2_score(ground_truth, predictions)
    metrics['NRMSE'] = metrics['RMSE'] / np.sqrt(range)
    metrics['MAPE'] = mean_absolute_percentage_error(ground_truth, predictions)
    return metrics

def linear_regression_fit(X_train, y_train):
  # print(len(features), "Training datapoints")
  # X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
  model = LinearRegression()
  # Fit the model to the training data
  model.fit(X_train, y_train)
  # Make predictions on the testing data
  return model

def linear_regression_predict(model, X_test):
  y_pred = model.predict(X_test)
  print("Feature Coefficients:")
  # print(X_train.columns)
  for feature, coef in zip(X_test.columns, model.coef_):
    print(f"{feature}: {coef}")
  return pd.DataFrame(y_pred, columns=['Predictions'])

def write_predictions(predictions, ground_truth=None, filename="predictions.csv"):
  # Create the results folder if it doesn't exist
  if not os.path.exists(results_folder):
    os.makedirs(results_folder)

  # Create the CSV file for the training predictions and ground truth
  predictions_file = os.path.join(results_folder, filename)

  # Write the training predictions and ground truth dfs to one csv file
  if ground_truth is not None:
    # concatenate the predictions and ground truth pandas dataframes
    predictions = pd.concat([predictions, ground_truth], axis=1)


  predictions.to_csv(predictions_file, index=False)

def write_metrics(metrics, filename):
  # Write the training error metrics and testing error metrics to a separate file
  metrics_file = os.path.join(results_folder, filename)
  with open(metrics_file, "w") as f:
    f.write("Test Error Metrics:\n")
    f.write("Mean Absolute Error (MAE): {}\n".format(metrics['MAE']))
    f.write("Mean Squared Error (MSE): {}\n".format(metrics['MSE']))
    f.write("Root Mean Squared Error (RMSE): {}\n".format(metrics['RMSE']))
    f.write("R-squared (R²) Score: {}\n".format(metrics['R2']))
    f.write("Normalized Root Mean Squared Error: {}\n".format(metrics['NRMSE']))
    f.write("Mean Absolute Percentage Error: {}\n".format(metrics['MAPE']))

def wind_to_power(df):
  return df * CAPACITY

def main():
#   filename = "3.6mw_wind_data.csv"
#   X, y = process_training_data(filename)
#   X_train, X_test, y_train, y_test = train_test_split(
#       X, y, test_size=0.2, random_state=0)
#   y_test = pd.DataFrame(y_test)
#   y_test.reset_index(drop=True, inplace=True)

#   model = random_forest.train_model(X_train, y_train)
#   y_pred = random_forest.predict(model, X_test)

#   # pdef = memory_problem_def(train_df, test_df)
#   # automl(pdef, test_df)

# #   model = linear_regression_fit(X_train, y_train)
# #   y_pred = linear_regression_predict(model, X_test)

# #   # model = lstm_fit(X_train, y_train)
# #   # y_pred = lstm_predict(model, X_test)

#   write_predictions(y_pred, y_test, "wind_predictions.csv")
#   test_metrics = evaluate_model(y_test, y_pred)
# #   print("Mean Absolute Error (MAE):", test_metrics['MAE'])
# #   print("Mean Squared Error (MSE):", test_metrics['MSE'])
# #   print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
# #   print("R-squared (R²) Score:", test_metrics['R2'])
# #   print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
# #   print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
#   write_metrics(test_metrics, 'wind_metrics_rf.txt')

#   # model.save(f"results/wind_model_rf.keras")
#   dump(model, "results/wind_model_rf.joblib")
#   # # load the model
  model = load("results/wind_model_rf.joblib")
#   model = tf.keras.models.load_model("results/wind_model_rf.keras")
  forecast = pull_weather_forecast("Pittsburgh, PA, US")
  forecast_df = format_wind_forecast(forecast)
  print(forecast_df.head())
  predictions = random_forest.predict(model, forecast_df)
#   # predictions = load_pipeline_and_predict(pdef, forecast_df)
  print(predictions)
# #   # print(type(forecast_df))

#   predictions = linear_regression_predict(model, forecast_df)

#   # predictions = lstm_predict(model, forecast_df)
#   print(predictions.head())


if __name__=="__main__":
    main()
