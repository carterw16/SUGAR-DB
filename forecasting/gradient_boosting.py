import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler, RobustScaler

def preprocess_data(load_df, scaler_x, scaler_y):
    # Ensure 'Time' is in datetime format and set as index
    load_df['Time'] = pd.to_datetime(load_df['Time'])
    load_df.set_index('Time', inplace=True)
    df = load_df

    # weather_df['Time'] = pd.to_datetime(weather_df['dt_iso'])
    # weather_df.set_index('Time', inplace=True)

    # print(load_df.head())
    # print(weather_df.head())
    # # Resample to get daily average load and temperature
    # daily_load = load_df['AverageLoad'].resample('D').mean()
    # daily_weather = weather_df.resample('D').mean()
    # df = pd.merge(load_df, weather_df, left_index=True, right_index=True, how='inner')

    # # Merge datasets on the daily index
    # df = pd.concat([daily_load, daily_weather], axis=1)

    # Drop any rows with NaN values that resulted from resampling
    df.dropna(inplace=True)

    # Extract temporal features from the index
    df['hour'] = df.index.hour
    df['day_of_week'] = df.index.dayofweek
    # df['month'] = df.index.month
    # df['day_of_month'] = df.index.day

    # Normalize weather features
    # scaler = RobustScaler()
    features = ['hour']
    # print(df.head())
    # df[features] = scaler_x.fit_transform(df[features].values)
    df['AverageLoad'] = df['AverageLoad'] / scaler_y
    df.reset_index(inplace=True)
    df.drop(columns=['Time'], inplace=True)

    # Separate features and target
    X = df.drop('AverageLoad', axis=1)
    y = df['AverageLoad']

    return X, y

# Train Model Function
def train_model(X_train, y_train):
    # Initialize the model with your chosen hyperparameters
    model = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
    # Train the model
    model.fit(X_train, y_train)
    return model

# Predict Model Function
def predict_model(model, X_new):
    predictions = model.predict(X_new)
    predictions_df = pd.DataFrame({
        # 'Hour': X_new['hour'],
        'Prediction': predictions
    })
    return predictions_df
