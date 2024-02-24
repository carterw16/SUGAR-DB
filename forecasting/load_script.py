import pandas as pd
import os
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_absolute_percentage_error
from wind_script import evaluate_model
import numpy as np
import requests
import urllib.parse
import time

def process_data(DATA_DIR):
  data_path = os.path.join(DATA_DIR, "demand_40.4396_-79.9839.csv")
  df = pd.read_csv(data_path)

  df = df.drop(columns=['time'])
  print(df.head())
  return df

def main():
  DATA_DIR = os.path.abspath('./data')
  df = process_data(DATA_DIR)

if __name__=="__main__":
    main()