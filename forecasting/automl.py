
import sys, os
import pandas as pd
PATH_TO_PICARD="/Users/carterweaver/Desktop/Summer2023/Auton/ngautonml"
sys.path.append(PATH_TO_PICARD)

from ngautonml.wrangler.wrangler import Wrangler
from ngautonml.wrangler.dataset import DatasetKeys, RoleName
from ngautonml.problem_def.problem_def  import ProblemDefinition
from ngautonml.algorithms.impl.algorithm_auto import FakeAlgorithmCatalogAuto
from ngautonml.wrangler.saver import Saver
from ngautonml.problem_def.output_config import OutputConfig

WIND_TARGET = "LV ActivePower (%)"

# MEMORY_PROBLEM_DEF = '''
# {
#     "dataset": {
#         "config": "memory",
#         "column_roles": {
#             "target": {
#                 "name": "target",
#             }
#         }
#     },
#     "problem_type": {
#         "task": "regression"
#     },
#     "cross_validation": {
#         "k": 10
#     },
#     "metrics" : {
#         "root_mean_squared_error" : {}
#     },
#     "hyperparams": ["disable_grid_search"]
# }
# '''

pdef_dict = {
  "dataset": {
    "config": "memory",
    "column_roles": {
      "target" : {
        "name": WIND_TARGET
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
  path = os.path.abspath('./automl_models/wind')
  # tmp_dir = path + '/wind'
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
  # print(problem_def.dataset_config.cols_with_role(role=RoleName.TARGET))
  te_ground_truth = wrangler.dataset(data=test_df[[WIND_TARGET]], key=DatasetKeys.GROUND_TRUTH)
  # te_ground_truth = wrangler.dataset(data=test_data.dataframe[[target]], key=DatasetKeys.GROUND_TRUTH)
  test_rankings = wrangler.rank(results=test_result, ground_truth=te_ground_truth)

  print(train_rankings)
  print(test_rankings)
  # save the best executable pipeline as a runnable model
  saver = Saver(output_config=output_config())
  saver.save_models(pipelines=test_result.executable_pipelines)
  # saver.save_all_models(test_result)
