from sklearn.preprocessing import MinMaxScaler
from keras.models import Sequential
from keras.layers import LSTM, Dense, Input
import pandas as pd
import numpy as np

def lstm_normalize(X):
  scaler = MinMaxScaler()
  X_scaled = scaler.fit_transform(X)
  X_lstm = X_scaled.reshape(
      (X_scaled.shape[0], 1, X_scaled.shape[1]))
  return X_lstm

def lstm_fit(X_train, y_train):
  print(len(X_train), "Training datapoints")
  # Split the dataset into training and testing sets
  # X_train, X_test, y_train, y_test = train_test_split(
  #     features, target, test_size=0.2, random_state=42)
  # Normalize the features
  # scaler = MinMaxScaler()
  X_train_lstm = lstm_normalize(X_train)
  # X_train_lstm = X_train
  # Train the LSTM model
  model = Sequential()
  model.add(Input(shape=(X_train_lstm.shape[1], X_train_lstm.shape[2])))
  model.add(LSTM(50))
  model.add(Dense(1))
  model.compile(loss='mean_squared_error', optimizer='adam')
  model.fit(X_train_lstm, y_train, epochs=50, batch_size=72, verbose=2)

  return model

def lstm_predict(model, input):
  input_lstm = lstm_normalize(input)
  y_pred = model.predict(input_lstm)
  # pull negative values to 0
  y_pred = np.maximum(y_pred, 0)
  return pd.DataFrame(y_pred, columns=['Predictions'])
