import pandas as pd
import os, datetime
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from keras.models import Sequential
from keras.layers import LSTM, Dense, Input
import requests
import config
import tensorflow as tf

WEATHER_API_KEY = config.WEATHER_API_KEY
DATA_DIR = os.path.abspath('./data')
results_folder = "results"

def process_training_data(filename):
  data_path = os.path.join(DATA_DIR, filename)
  full_df = pd.read_csv(data_path)

  cap = 3600
  full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / cap
  full_df['Time'] = pd.to_datetime(full_df['Date/Time'], format='%d %m %Y %H:%M')
  # load_df = load_df[['Time','Aggregate']]
  # full_df = full_df.dropna(axis=0, how='any')
  # drop the Date/time column from full_df
  full_df.drop(columns=['Date/Time'], inplace=True)
  full_df.set_index('Time', inplace=True)
  full_df = full_df.apply(pd.to_numeric)
  full_df = full_df.resample('h').mean().fillna(0)
  full_df.reset_index(inplace=True)

  # weather_df = pd.read_csv(weather_path)
  # weather_df = weather_df[['dt', 'temp', 'pressure', 'humidity', 'clouds_all']]
  # weather_df['dt'] = pd.to_datetime(weather_df['dt'], unit='s')
  # weather_df=weather_df.rename(columns={'dt': 'Time'})

  # df = full_df.merge(weather_df, on='Time')

  full_df['month'] = full_df['Time'].dt.month
  full_df['day'] = full_df['Time'].dt.dayofweek
  full_df['hour'] = full_df['Time'].dt.hour
  X = full_df[["month","day","hour","Wind Speed (m/s)"]]
  y = full_df[["LV ActivePower (%)"]]
  # print(full_df.head())
  # print(X.head())
  # print(y.head())
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
      # convert timestamp object to datetime object
      # timestamp = timestamp.to_pydatetime()
      # convert GMT timestamp to local time for forecastes location
      # timestamp = timestamp + datetime.timedelta(hours=forecast["timezone"]/3600)
      # new Date((data.sys.sunrise + data.timezone) * 1000)

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

def linear_regression(X_train, X_test, y_train, y_test):
  # print(len(features), "Training datapoints")
  # X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)
  model = LinearRegression()
  # Fit the model to the training data
  model.fit(X_train, y_train)
  # Make predictions on the testing data
  y_pred = model.predict(X_test)
  print("Feature Coefficients:")
  print(X_train.columns)
  for feature, coef in zip(X_train.columns, model.coef_):
    print(f"{feature}: {coef}")
  # Evaluate the model
  metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", metrics['MAE'])
  print("Mean Squared Error (MSE):", metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", metrics['RMSE'])
  print("R-squared (R²) Score:", metrics['R2'])
  print("Normalized Root Mean Squared Error:", metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", metrics['MAPE'])
  print(y_pred, y_test)
  return metrics

def lstm_normalize(X):
  scaler = MinMaxScaler()
  X_scaled = scaler.fit_transform(X)
  X_lstm = X_scaled.reshape(
      (X_scaled.shape[0], 1, X_scaled.shape[1]))
  return X_lstm

def lstm_fit(X_train, y_train):
  print(len(X_train), "Training datapoints")
  # Split the dataset into training and testing sets
  # X_train, X_test, y_train, y_test = train_test_split(
  #     features, target, test_size=0.2, random_state=42)
  # Normalize the features
  # scaler = MinMaxScaler()
  X_train_lstm = lstm_normalize(X_train)
  # Train the LSTM model
  model = Sequential()
  model.add(Input(shape=(X_train_lstm.shape[1], X_train_lstm.shape[2])))
  model.add(LSTM(50))
  model.add(Dense(1))
  model.compile(loss='mean_squared_error', optimizer='adam')
  model.fit(X_train_lstm, y_train, epochs=50, batch_size=72, verbose=2)

  return model

def lstm_predict(model, input):
  input_lstm = lstm_normalize(input)
  y_pred = model.predict(input_lstm)
  # return one column dataframe containg the predictions
  return pd.DataFrame(y_pred, columns=['Predictions'])

  # return pd.DataFrame(y_pred)

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

def main():
  filename = "3.6mw_wind_data.csv"
  X, y = process_training_data(filename)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42,  shuffle=False)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  # linear_regression(X_train, X_test, y_train, y_test)
  model = lstm_fit(X_train, y_train)
  # Make predictions on the testing data
  y_pred = lstm_predict(model, X_test)
  write_predictions(y_pred, y_test, "wind_predictions.csv")

  test_metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", test_metrics['MAE'])
  print("Mean Squared Error (MSE):", test_metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
  print("R-squared (R²) Score:", test_metrics['R2'])
  print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
  write_metrics(test_metrics, 'wind_metrics_lstm.txt')

  model.save(f"results/wind_model_lstm.keras")
  # load the model
  model = tf.keras.models.load_model("results/wind_model_lstm.keras")
  forecast = pull_weather_forecast("Pittsburgh, PA, US")
  forecast_df = format_wind_forecast(forecast)
  print(forecast_df.head())
  # print(type(forecast_df))

  predictions = lstm_predict(model, forecast_df)
  print(predictions.head())


if __name__=="__main__":
    main()