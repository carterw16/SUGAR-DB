import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import linear_regression, lstm_fit, lstm_predict, write_metrics, write_predictions, evaluate_model
import numpy as np
import requests
import urllib.parse
import time
import config
# SOLAR_API_KEY = os.getenv("SOLAR_API_KEY")
# EMAIL = "carterw@andrew.cmu.edu"
# BASE_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-5min-download.json?"
# POINTS = [
# '2539138'
# ]
WEATHER_API_KEY = config.WEATHER_API_KEY
DATA_DIR = os.path.abspath('./data')
results_folder = "results"
# STC solar irradiation for PV modules in W/m^2
GHI_STANDARD = 1000

def process_training_data(filename):
  data_path = os.path.join(DATA_DIR, filename)
  df = pd.read_csv(data_path)
  df = df.dropna(axis=1, how='any')
  df = df.astype('float64')
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

def main():
  filename = 'solarrad_40.44_-79.99_2022.csv'
  X, y = process_training_data(filename)
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, test_size=0.2, random_state=42)
  y_test = pd.DataFrame(y_test)
  y_test.reset_index(drop=True, inplace=True)
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")

  # linear_regression(X,y)
  model = lstm_fit(X_train, y_train)
  # Make predictions on the testing data
  y_pred = lstm_predict(model, X_test)
  write_predictions(y_pred, y_test, "solar_predictions.csv")

  test_metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", test_metrics['MAE'])
  print("Mean Squared Error (MSE):", test_metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", test_metrics['RMSE'])
  print("R-squared (R²) Score:", test_metrics['R2'])
  print("Normalized Root Mean Squared Error:", test_metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", test_metrics['MAPE'])
  write_metrics(test_metrics, 'solar_metrics_lstm.txt')

  model.save(f"results/solar_model_lstm.keras")

if __name__=="__main__":
    main()