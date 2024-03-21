import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from keras.models import Sequential
from keras.layers import LSTM, Dense, Input


DATA_DIR = os.path.abspath('./data')

def process_data(filename):
  data_path = os.path.join(DATA_DIR, filename)
  full_df = pd.read_csv(data_path)
  cap = 3600
  full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / cap
  X = full_df[["Wind Speed (m/s)","Wind Direction (°)"]]
  y = full_df["LV ActivePower (%)"]
  return X, y

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

def linear_regression(features, target):
  print(len(features), "Training datapoints")
  # Create baseline model with linear regression
  # Split the dataset into training and testing sets
  X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)

  # Initialize the linear regression model
  model = LinearRegression()

  # Fit the model to the training data
  model.fit(X_train, y_train)

  # Make predictions on the testing data
  y_pred = model.predict(X_test)

  # Print the coefficients of the linear regression model
  # print("Coefficients:", model.coef_)
  # Show which coefficient is associated with which feature
  print("Feature Coefficients:")
  for feature, coef in zip(features.columns, model.coef_):
    print(f"{feature}: {coef}")

  # Evaluate the model
  metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", metrics['MAE'])
  print("Mean Squared Error (MSE):", metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", metrics['RMSE'])
  print("R-squared (R²) Score:", metrics['R2'])
  print("Normalized Root Mean Squared Error:", metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", metrics['MAPE'])

  return metrics

def lstm_regression(features, target):
  print(len(features), "Training datapoints")
  # Split the dataset into training and testing sets
  X_train, X_test, y_train, y_test = train_test_split(
      features, target, test_size=0.2, random_state=42)
  # Step 5: Normalize the features
  scaler = MinMaxScaler()
  X_train_scaled = scaler.fit_transform(X_train)
  X_test_scaled = scaler.transform(X_test)

  # Step 6: Train the LSTM model
  X_train_lstm = X_train_scaled.reshape((X_train_scaled.shape[0], 1, X_train_scaled.shape[1]))
  model = Sequential()
  model.add(Input(shape=(X_train_lstm.shape[1], X_train_lstm.shape[2])))
  model.add(LSTM(50))
  model.add(Dense(1))
  model.compile(loss='mean_squared_error', optimizer='adam')
  model.fit(X_train_lstm, y_train, epochs=50, batch_size=72, verbose=2)

  # Step 7: Make predictions
  X_test_lstm = X_test_scaled.reshape((X_test_scaled.shape[0], 1, X_test_scaled.shape[1]))
  y_pred = model.predict(X_test_lstm)

  # Print the first 5 predictions
  print(y_pred[:5])
  print(y_test[:5])

  # Step 8: Evaluate the model
  metrics = evaluate_model(y_test, y_pred)
  print("Mean Absolute Error (MAE):", metrics['MAE'])
  print("Mean Squared Error (MSE):", metrics['MSE'])
  print("Root Mean Squared Error (RMSE):", metrics['RMSE'])
  print("R-squared (R²) Score:", metrics['R2'])
  print("Normalized Root Mean Squared Error:", metrics['NRMSE'])
  print("Mean Absolute Percentage Error:", metrics['MAPE'])

  return metrics

def main():
  filename = "3.6mw_wind_data.csv"
  X, y = process_data(filename)
  linear_regression(X,y)
  # lstm_regression(X,y)

if __name__=="__main__":
    main()