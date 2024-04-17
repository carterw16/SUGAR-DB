import sys, os
import pandas as pd
from pathlib import Path
PATH_TO_PICARD="/Users/carterweaver/Desktop/Summer2023/Auton/ngautonml"
sys.path.append(PATH_TO_PICARD)

from ngautonml.wrangler.wrangler import Wrangler
from ngautonml.wrangler.dataset import DatasetKeys, RoleName
from ngautonml.problem_def.problem_def  import ProblemDefinition
from ngautonml.algorithms.impl.algorithm_auto import AlgorithmCatalogAuto
from ngautonml.wrangler.saver import Saver
from ngautonml.problem_def.output_config import OutputConfig
from ngautonml.instantiator.instantiator_factory import InstantiatorFactory
from ngautonml.instantiator.json_loader import JsonLoader
from ngautonml.instantiator.json_instantiator import JsonInstantiator
from ngautonml.executor.executor_kind import ExecutorKind
from load_script import process_training_data, format_load_forecast
from wind_script import pull_weather_forecast

WIND_TARGET = "LV ActivePower (%)"
LOAD_TARGET = "AverageLoad"

pdef_dict = {
  "dataset": {
    "config": "memory",
    "column_roles": {
      "target" : {
        "name": LOAD_TARGET
      }
    }
  },
  "problem_type": {
    "task": "regression"
  },
  "cross_validation": {
    "k": 5
  },
  "output": {},
  "metrics": {
    "r2_score" : {},
    "root_mean_squared_error": {}
  },
  "hyperparams": ["disable_grid_search"]
}

def memory_problem_def(train_df: pd.DataFrame,
                      test_df: pd.DataFrame) -> ProblemDefinition:
  problem_def = pdef_dict

  return ProblemDefinition(
      problem_def=problem_def,
      train_df=train_df,
      test_df=test_df)

def output_config() -> OutputConfig:
  path = os.path.abspath('./automl_models/load')
  # tmp_dir.mkdir()
  if not os.path.exists(path):
    os.makedirs(path)
  config = OutputConfig(clause={
      'path': path,
      'instantiations': [
          'stub_executor_kind',
          'simple'
      ],
      'file_type': 'csv'
  })
  return config

def automl(problem_def, test_df):
  wrangler = Wrangler(problem_definition=problem_def)
  wrangler_result = wrangler.fit_predict_rank()

  train_result = wrangler_result.train_results
  train_rankings = wrangler_result.rankings
  test_result = wrangler_result.test_results
  te_ground_truth = wrangler.dataset(data=test_df[[LOAD_TARGET]], key=DatasetKeys.GROUND_TRUTH)
  test_rankings = wrangler.rank(results=test_result, ground_truth=te_ground_truth)
  best_test_pipeline = test_rankings["r2_score"].best(1)[0].pipeline_des
  bound_pipeline = test_result.bound_pipelines[best_test_pipeline]
  print(train_rankings["r2_score"])
  print(test_rankings["r2_score"])
  # print(te_ground_truth)
  # print(test_result.predictions)
  # save the best executable pipeline as a runnable model
  saver = Saver(output_config=output_config())
  model_paths = saver.save_models(pipelines=test_result.executable_pipelines)
  json_saver = JsonInstantiator(saver=saver)
  json_saver.save(bound_pipeline, model_paths=model_paths)
  # saver.save_models(pipelines=test_result.executable_pipelines)
  # saver.save_all_models(test_result)

def load_pipeline_and_predict(pdef, input_data, pipeline_file="./automl_models/wind/pipelines/tabular_regression@sklearn.ensemble.randomforestregressor:max_depth=none,min_samples_split=2.json"):
  wrangler = Wrangler(problem_definition=pdef)
  loader = JsonLoader(
        saver_version='1.0',
        algorithm_catalog=AlgorithmCatalogAuto(),
        pipeline_file=Path(pipeline_file),
        load_models=True)
  bound_pipeline = loader.pipeline
  instantiator = InstantiatorFactory().build(kind=ExecutorKind('simple'))
  pipelines = {bound_pipeline.designator: instantiator.instantiate(bound_pipeline)}
  test_data = input_data

  test_dataset = wrangler.dataset(data=test_data)

  predictions = wrangler.predict(new_data=test_dataset, trained_pipelines=pipelines)
  return predictions

def main():
  load_file = 'average_loads.csv'
  weather_file = 'openweather_hourly_2013_2023/loughborough.csv'
  train_df, test_df = process_training_data(load_file, weather_file)
  print(train_df)
  print(test_df)
  pdef = memory_problem_def(train_df, test_df)
  automl(pdef, test_df)

  # forecast = pull_weather_forecast("Pittsburgh, PA, US")
  # load_df = format_load_forecast(forecast)
  # print(load_df)
  # pipeline_file = "./automl_models/load/pipelines/tabular_regression@sklearn.ensemble.randomforestregressor:max_depth=none,min_samples_split=2.json"
  # predictions = load_pipeline_and_predict(pdef, test_df, pipeline_file)
  # print(predictions)

if __name__ == "__main__":
  main()