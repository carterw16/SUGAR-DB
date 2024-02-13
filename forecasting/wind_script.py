import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
import numpy as np

def process_data(DATA_DIR):
  data_path = os.path.join(DATA_DIR, "3.6mw_wind_data.csv")
  full_df = pd.read_csv(data_path)
  cap = 3600
  full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / cap
  # train_df, test_df = train_test_split(full_df, test_size=0.2)
  return full_df

def evaluate_model(ground_truth, predictions):
    variance = np.var(ground_truth)
    metrics = {}
    metrics['MAE'] = mean_absolute_error(ground_truth, predictions)
    metrics['MSE'] = mean_squared_error(ground_truth, predictions)
    metrics['RMSE'] = np.sqrt(metrics['MSE'])
    metrics['R2'] = r2_score(ground_truth, predictions)
    metrics['RELMSE'] = metrics['RMSE'] / np.sqrt(variance)
    metrics['MAPE'] = mean_absolute_percentage_error(ground_truth, predictions)
    return metrics

def baseline(df):
  target = df["LV ActivePower (%)"]
  # Create baseline model with persistence algorithm
  predictions = target.shift(1)
  predictions.iloc[0] = target[0]

  metrics = evaluate_model(target, predictions)
  # print("Mean Absolute Error (MAE):", mae)
  # print("Mean Squared Error (MSE):", mse)
  # print("Root Mean Squared Error (RMSE):", rmse)
  # print("R-squared (R²) Score:", r2)
  print("Relative Mean Squared Error:", metrics['RELMSE'])
  print("Mean Absolute Percentage Error:", metrics['MAPE'])

  return metrics

def main():
  DATA_DIR = os.path.abspath('./data')
  df = process_data(DATA_DIR)
  baseline(df)

if __name__=="__main__":
    main()