import os
import sys
# parent_dir = os.getcwd()
# print(parent_dir)
sys.path.append(os.getcwd() + "/forecasting")
sys.path.append(os.getcwd() + "/sugar")
from wind_script import *
from load_script import *
from solar_script import *
import config

from runSUGAR3 import *

def wind_forecast(location="Pittsburgh, PA, US"):
  # load the model
  wind_model = tf.keras.models.load_model("forecasting/results/wind_model_lstm.keras")
  forecast = pull_weather_forecast(location)
  wind_forecast_df = format_wind_forecast(forecast)
  # print(wind_forecast_df.head())
  wind_predictions = lstm_predict(wind_model, wind_forecast_df)
  # add hour column to predictions
  wind_predictions['hour'] = wind_forecast_df.hour
  # print(wind_predictions)
  return wind_predictions

def solar_forecast(location="Pittsburgh, PA, US"):
  # load the model
  solar_model = tf.keras.models.load_model("forecasting/results/solar_model_lstm.keras")
  forecast = pull_weather_forecast(location)
  solar_forecast_df = format_solar_forecast(forecast)
  # print(solar_forecast_df.head())
  solar_predictions = lstm_predict(solar_model, solar_forecast_df)
  # print(solar_predictions)
  solar_predictions['hour'] = solar_forecast_df.Hour
  # shift all predictions back 3 hours but still start at the same hour
  solar_predictions['hour'] = (solar_predictions['hour'] + 21) % 24
  # remove first 3 hours
  solar_predictions = solar_predictions.iloc[3:].reset_index(drop=True)
  # start
  return solar_predictions

def load_forecast(location="Pittsburgh, PA, US"):
  # load the model
  load_model = tf.keras.models.load_model("forecasting/results/load_model_lstm.keras")
  forecast = pull_weather_forecast(location)
  load_forecast_df = format_load_forecast(forecast)
  # print(load_forecast_df.head())
  load_predictions = lstm_predict(load_model, load_forecast_df)
  # print(load_predictions)
  load_predictions['hour'] = load_forecast_df.hour
  return load_predictions

def run_all():
  location = "Austin, TX, US"
  #predictions = solar_forecast(location)
  data = {
  "Predictions": np.linspace(0,1,24),
  "hour": list(range(1, 25))
  }

  predictions = pd.DataFrame(data)
  print(predictions)

  # run SUGAR3
  settings['multi settings']['PV capacity'] = predictions['Predictions'][0:PERIODS]
  settings['multi settings']['wind capacity'] = np.linspace(1,1,PERIODS)
  main(case, settings, features)

  


if __name__=="__main__":
  run_all()
