import matplotlib.pyplot as plt
import seaborn as sns

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