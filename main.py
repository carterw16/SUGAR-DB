import os
import sys
# parent_dir = os.getcwd()
# print(parent_dir)
from pathlib import Path
parent_dir = str(Path(__file__).parent)
sys.path.append(parent_dir + "/forecasting")
sys.path.append(parent_dir + "/sugar")
# from wind_script import *
# from solar_script import *
import config

from joblib import dump, load
from wind_script import pull_weather_forecast, format_wind_forecast
from load_script import format_load_forecast
from solar_script import format_solar_forecast, df_ghi_to_power
import random_forest

# from sugar.runSUGAR3 import *


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
  model = load(f"{parent_dir}/forecasting/results/solar_model_rf.joblib")
  forecast = pull_weather_forecast(location)
  solar_forecast_df, solar_intensity = format_solar_forecast(forecast)
  solar_predictions = random_forest.predict(model, solar_forecast_df)
  # convert ghi to power
  power_predictions = df_ghi_to_power(solar_predictions)
  # power_predictions = solar_predictions.apply(lambda x: ghi_to_power_factor(x*GHI_STANDARD, SOLAR_CAPACITY, T_COEFF, T_STANDARD))
  # power_predictions = solar_predictions
  power_predictions['hour'] = solar_forecast_df.Hour
  # shift all predictions back 3 hours but still start at the same hour
  power_predictions['hour'] = (power_predictions['hour'] + 21) % 24
  # remove first 3 hours
  power_predictions = power_predictions.iloc[3:].reset_index(drop=True)
  power_predictions = power_predictions[:24]
  # print(power_predictions)
  print(solar_intensity)
  power_predictions["Predictions"] *= solar_intensity["SI"]
  return power_predictions

def load_forecast(location="Pittsburgh, PA, US"):
  load_model = load(f"{parent_dir}/forecasting/results/load_model_rf.joblib")

  forecast = pull_weather_forecast(location)
  load_forecast_df = format_load_forecast(forecast)
  load_predictions = random_forest.predict(load_model, load_forecast_df)
  load_predictions['hour'] = load_forecast_df.hour
  load_predictions = load_predictions[:24]
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
  solar_pred = solar_forecast("Pittsburgh, PA, US")
  wind_pred = wind_forecast("Pittsburgh, PA, US")
  print(wind_pred)
  print(solar_pred)
  print(prediction)
