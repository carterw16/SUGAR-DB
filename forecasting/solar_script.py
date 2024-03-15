import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import evaluate_model, linear_regression
import numpy as np
import requests
import urllib.parse
import time
from dotenv import load_dotenv
load_dotenv()
SOLAR_API_KEY = os.getenv("SOLAR_API_KEY")
EMAIL = "carterw@andrew.cmu.edu"
BASE_URL = "https://developer.nrel.gov/api/nsrdb/v2/solar/psm3-5min-download.json?"
POINTS = [
'2539138'
]

def process_data(DATA_DIR):
  data_path = os.path.join(DATA_DIR, "solarrad_40.44_-79.99_2022.csv")
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

def main():
  DATA_DIR = os.path.abspath('./data')
  df = process_data(DATA_DIR)
  X = df.drop(columns=['Year','GHI'])
  y = df['GHI']
  linear_regression(X, y)

if __name__=="__main__":
    main()