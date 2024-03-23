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
from wind_script import evaluate_model, linear_regression, lstm_fit, lstm_predict, write_metrics, write_predictions
import config
WEATHER_API_KEY = config.WEATHER_API_KEY
DATA_DIR = os.path.abspath('./data')

def process_training_data(load_file, weather_file):
  load_path = os.path.join(DATA_DIR, load_file)
  weather_path = os.path.join(DATA_DIR, weather_file)
  load_df = pd.read_csv(load_path)

  load_df['Time'] = pd.to_datetime(load_df['Time'])
  load_df = load_df[['Time','Aggregate']]
  load_df = load_df.dropna(axis=0, how='any')
  load_df['Time'] = pd.to_datetime(load_df['Time'])
  load_df.set_index('Time', inplace=True)
  load_df = load_df.resample('h').mean().fillna(0)
  load_df.reset_index(inplace=True)

  weather_df = pd.read_csv(weather_path)
  weather_df = weather_df[['dt', 'temp', 'pressure', 'humidity', 'clouds_all']]
  weather_df['dt'] = pd.to_datetime(weather_df['dt'], unit='s')
  weather_df=weather_df.rename(columns={'dt': 'Time'})

  df = load_df.merge(weather_df, on='Time')

  df['month'] = df['Time'].dt.month
  df['day'] = df['Time'].dt.dayofweek
  df['hour'] = df['Time'].dt.hour
  X = df.drop(columns=['Aggregate', 'Time'])
  y = df['Aggregate']
  return X,y

# def download_weather_data():
#   # Define the location and time range
#   city = "Loughborough"
#   country_code = "GB"
#   start_date = "2013-10-09"
#   end_date = "2019-12-31"
#   latitude = 40.4406  # Latitude of Pittsburgh
#   longitude = -79.9959  # Longitude of Pittsburgh
#   cnt = 24  # Number of hourly data points to retrieve per day

#   # Convert start_date to Unix timestamp (seconds since Epoch)
#   start_timestamp = int(datetime.strptime(start_date, "%Y-%m-%d").timestamp())

#   # Make the API request
#   base_url = "https://history.openweathermap.org/data/2.5/history/city"
#   url = f"{base_url}?q={city},{country_code}&type=hour&start={start_timestamp}&cnt={cnt}&appid={WEATHER_API_KEY}"

#   # Send the request
#   response = requests.get(url)

#   # Check if the request was successful (status code 200)
#   if response.status_code == 200:
#     # Extract the hourly weather data from the response
#     weather_data = response.json()
#     print("success")
#     # Print the hourly weather data
#     for data_point in weather_data["list"]:
#       timestamp = datetime.utcfromtimestamp(data_point["dt"]).strftime('%Y-%m-%d %H:%M:%S')
#       temperature = data_point["main"]["temp"]
#       humidity = data_point["main"]["humidity"]
#       wind_speed = data_point["wind"]["speed"]
#       print(f"Timestamp: {timestamp}, Temperature: {temperature}°C, Humidity: {humidity}%, Wind Speed: {wind_speed} m/s")

#   else:
#     print(f"Error: {response.status_code}, {response.text}")

def main():
  load_file = 'UK_house_loads/House_1.csv'
  weather_file = 'openweather_hourly_2013_2023/loughborough.csv'
  X,y = process_training_data(load_file, weather_file)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")

  # linear_regression(X,y)
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

if __name__=="__main__":
  main()