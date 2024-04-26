import os
import sys
# parent_dir = os.getcwd()
# print(parent_dir)
from pathlib import Path
parent_dir = os.getcwd()

sys.path.append(os.getcwd() + "/forecasting")
sys.path.append(os.getcwd() + "/sugar")
from wind_script import *
from load_script import *
from solar_script import *
import config

# from sugar.runSUGAR3 import *

GHI_STANDARD = 1000
T_STANDARD = 25
SOLAR_CAPACITY = 400
T_COEFF = -0.003

def wind_forecast(location="Pittsburgh, PA, US"):
  # load the model
  wind_model = load(f"{parent_dir}/forecasting/results/wind_model_rf.joblib")
  forecast = pull_weather_forecast(location)
  wind_forecast_df = format_wind_forecast(forecast)
  wind_predictions = random_forest.predict(wind_model, wind_forecast_df)
  wind_predictions['hour'] = wind_forecast_df.hour
  wind_predictions = wind_predictions[:24]
  return wind_predictions

def solar_forecast(location="Pittsburgh, PA, US"):
  # load the model
  # solar_model = tf.keras.models.load_model("forecasting/results/solar_model_lstm.keras")
  model = load("forecasting/results/solar_model_rf.joblib")
  # forecast = pull_weather_forecast("Pittsburgh, PA, US")
  # forecast_df = format_solar_forecast(forecast)
  # print(forecast_df.head())
  # print(type(forecast_df))

  # predictions = lstm_predict(model, forecast_df)

  forecast = pull_weather_forecast(location)
  solar_forecast_df = format_solar_forecast(forecast)
  # print(solar_forecast_df.head())
  # solar_predictions = lstm_predict(solar_model, solar_forecast_df)
  solar_predictions = random_forest.predict(model, solar_forecast_df)
  # power_predictions = solar_predictions.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, SOLAR_CAPACITY, T_COEFF, T_STANDARD))
  power_predictions = solar_predictions
  # print(solar_predictions)
  power_predictions['hour'] = solar_forecast_df.Hour
  # shift all predictions back 3 hours but still start at the same hour
  power_predictions['hour'] = (power_predictions['hour'] + 21) % 24
  # remove first 3 hours
  power_predictions = power_predictions.iloc[3:].reset_index(drop=True)
  # start
  power_predictions = power_predictions[:24]
  return power_predictions

def load_forecast(location="Pittsburgh, PA, US"):
  load_model = load("forecasting/results/load_model_rf.joblib")

  # load_predictions = linear_regression_predict(model, load_df)
  # load_predictions = gradient_boosting.predict_model(model, load_df)

  # load the model
  forecast = pull_weather_forecast(location)
  load_forecast_df = format_load_forecast(forecast)
  load_predictions = random_forest.predict(load_model, load_forecast_df)
  load_predictions['hour'] = load_forecast_df.hour
  return load_predictions

def run_all():
  location = "Pittsburgh, PA, US"
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
  prediction = load_forecast("Pittsburgh, PA, US")
  print(prediction)
