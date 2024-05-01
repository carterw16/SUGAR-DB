import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler, RobustScaler

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
