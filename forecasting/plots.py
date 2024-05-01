import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.interpolate import make_interp_spline


def plot_heatmap(df):
  # Create a pivot table
  pivot_table = df.pivot_table(values='AverageLoad', index=df.index.hour, columns=df.index.month, aggfunc='mean')

  plt.figure(figsize=(12, 8))
  sns.heatmap(pivot_table, annot=True, fmt=".1f", linewidths=.5, cmap='Blues')
  plt.title('Average Load by Hour and Month')
  plt.xlabel('Month')
  plt.ylabel('Hour of Day')
  plt.show()

def plot_load_over_time(df):
  plt.figure(figsize=(15, 7))
  plt.plot(df['month'], df['AverageLoad'], label='Average Load')
  plt.title('Average Load of 20 UK Houses Over Time')
  plt.xlabel('Month')
  plt.ylabel('Average Load')
  plt.legend()
  plt.show()

def plot_load_vs_temp(df):
  # Calculate average daily load and average daily temperature
  daily_data = df.resample('D').agg({'AverageLoad':'mean', 'temp':'mean'})

  # Now create the scatter plot
  plt.figure(figsize=(10, 6))
  plt.scatter(daily_data['temp'], daily_data['AverageLoad'], alpha=0.5)
  plt.title('Load vs Average Daily Temperature')
  plt.xlabel('Average Daily Temperature (°C)')
  plt.ylabel('Average Daily Load')
  plt.grid(True)

  # Optionally, you can calculate and display the correlation coefficient
  correlation = daily_data['AverageLoad'].corr(daily_data['temp'])
  print(f"The correlation coefficient between daily load and average daily temperature is: {correlation:.2f}")
  plt.show()

def plot_hourly_load(df):
  hourly_load = df.groupby(df.index.hour).mean()

  plt.figure(figsize=(12, 6))
  hourly_load['AverageLoad'].plot(marker='o')
  plt.title('Average Hourly Load')
  plt.xlabel('Hour of the Day')
  plt.ylabel('Average Load')
  plt.xticks(range(0, 24))
  plt.grid(True)
  plt.show()

def plot_actual_vs_pred(dfs, labels):
  plt.figure(figsize=(11, 9))
  for df, label in zip(dfs, labels):
    # Interpolating for a smoother line
    xnew = np.linspace(df['Hour'].min(), df['Hour'].max(), 200)  # 300 represents number of points to make between min and max
    spl_actual = make_interp_spline(df['Hour'], df['Actual'], k=3)  # BSpline object
    smooth_actual = spl_actual(xnew)
    spl_predicted = make_interp_spline(df['Hour'], df['Predicted'], k=3)
    smooth_predicted = spl_predicted(xnew)
    plt.plot(xnew, smooth_actual, label='Actual Load on ' + label, linewidth=3)
    plt.plot(xnew, smooth_predicted, label='Predicted Load on ' + label, linewidth=3)

  plt.title('Load Predictions vs Actual Load Over 24 Hours - Test Set', fontsize=24)
  plt.ylim(0, .8)
  plt.xlim(0, 23)
  plt.xlabel('Hour of Day', fontsize=20)
  plt.ylabel('Load', fontsize=20)
  plt.legend(fontsize=18)
  plt.tight_layout()
  plt.show()