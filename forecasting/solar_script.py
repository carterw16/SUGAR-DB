import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
import numpy as np
import requests
import urllib.parse
import time

API_KEY = "HCbzILUfWXbR8c6x9YgjSIb6dUWkvGcQL0Yj2Gm2"
EMAIL = "carterw@andrew.cmu.edu"
BASE_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-5min-download.json?"
POINTS = [
'2539138'
]

def process_data(DATA_DIR):
  data_path = os.path.join(DATA_DIR, "2539138_40.44_-79.99_2022.csv")
  df = pd.read_csv(data_path)

  df.columns = df.iloc[1]
  df = df[2:]
  df = df.dropna(axis=1, how='any')
  df = df.astype('float64')
  # print(df.iloc[-1])
  # cap = 3600
  # full_df['LV ActivePower (%)'] = full_df['LV ActivePower (kW)'] / cap
  # train_df, test_df = train_test_split(full_df, test_size=0.2)
  return df

def evaluate_model(ground_truth, predictions):
  print(ground_truth.iloc[0].dtype)
  variance = np.var(ground_truth)
  metrics = {}
  metrics['MAE'] = mean_absolute_error(ground_truth, predictions)
  metrics['MSE'] = mean_squared_error(ground_truth, predictions)
  metrics['RMSE'] = np.sqrt(metrics['MSE'])
  metrics['R2'] = r2_score(ground_truth, predictions)
  metrics['RELMSE'] = metrics['RMSE'] / np.sqrt(variance)
  metrics['MAPE'] = mean_absolute_percentage_error(ground_truth, predictions)
  return metrics

def linear_regression(df):
  # Create baseline model with linear regression

  # Split the dataset into features (X) and target variable (y)
  X = df.drop(columns=['Year','GHI'])
  y = df['GHI']

  # Split the dataset into training and testing sets
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

  # Initialize the linear regression model
  model = LinearRegression()

  # Fit the model to the training data
  model.fit(X_train, y_train)

  # Make predictions on the testing data
  y_pred = model.predict(X_test)

  # Print the coefficients of the linear regression model
  print("Coefficients:", model.coef_)
  metrics = evaluate_model(y_test, y_pred)
  # print("Mean Absolute Error (MAE):", mae)
  print("Mean Squared Error (MSE):", metrics['MSE'])
  # print("Root Mean Squared Error (RMSE):", rmse)
  print("R-squared (R²) Score:", metrics['R2'])
  print("Relative Mean Squared Error:", metrics['RELMSE'])
  print("Mean Absolute Percentage Error:", metrics['MAPE'])

  return metrics

def download_solar_data():
  for year in ['2022']:
    print(f"Processing name: {year}")
    for id, location_ids in enumerate(POINTS):
      print(f'Making request for point group {id + 1} of {len(POINTS)}...')
      input_data = {
        'years': [year],
        'location_ids': location_ids,
        'attributes': 'dhi,dni,wind_speed,relative_humidity,air_temperature,surface_pressure,ghi,clearsky_ghi,clearsky_dni,clearsky_dhi,cloud_type',
        'interval': '30',
        'api_key': API_KEY,
        'email': EMAIL,
      }
      if '.csv' in BASE_URL:
        url = BASE_URL + urllib.parse.urlencode(data, True)
        # Note: CSV format is only supported for single point requests
        # Suggest that you might append to a larger data frame
        data = pd.read_csv(url)
        print(f'Response data (you should replace this print statement with your processing): {data}')
        # You can use the following code to write it to a file
        # data.to_csv('SingleBigDataPoint.csv')
      else:
        headers = {
          'x-api-key': API_KEY
        }
        data = get_response_json_and_handle_errors(requests.post(BASE_URL, input_data, headers=headers))
        download_url = data['outputs']['downloadUrl']
        # You can do with what you will the download url
        print(data['outputs']['message'])
        print(f"Data can be downloaded from this url when ready: {download_url}")

        # Delay for 1 second to prevent rate limiting
        time.sleep(1)
      print(f'Processed')

def get_response_json_and_handle_errors(response: requests.Response) -> dict:
  """Takes the given response and handles any errors, along with providing
  the resulting json
  Parameters
  ----------
  response : requests.Response
      The response object

  Returns
  -------
  dict
      The resulting json
  """
  if response.status_code != 200:
    print(f"An error has occurred with the server or the request. The request response code/status: {response.status_code} {response.reason}")
    print(f"The response body: {response.text}")
    exit(1)
  try:
    response_json = response.json()
  except:
    print(f"The response couldn't be parsed as JSON, likely an issue with the server, here is the text: {response.text}")
    exit(1)
  if len(response_json['errors']) > 0:
    errors = '\n'.join(response_json['errors'])
    print(f"The request errored out, here are the errors: {errors}")
    exit(1)
  return response_json

def main():
  DATA_DIR = os.path.abspath('./data')
  df = process_data(DATA_DIR)
  linear_regression(df)

if __name__=="__main__":
    main()