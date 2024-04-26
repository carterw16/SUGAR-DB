import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer

def train_model(X, y):
    """
    Trains the RandomForestRegressor model pipeline.

    Parameters:
    - X: DataFrame containing the features.
    - y: Series or array containing the target.

    Returns:
    - The trained model pipeline.
    """
    # Identifying categorical and numerical columns
    categorical_cols = X.select_dtypes(include=['object', 'category']).columns
    numerical_cols = X.select_dtypes(exclude=['object', 'category']).columns

    # Preprocessing for numerical data
    numerical_transformer = SimpleImputer(strategy='mean')

    # Preprocessing for categorical data
    categorical_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='most_frequent')),
        ('onehot', OneHotEncoder(handle_unknown='ignore', sparse_output=False))
    ])

    # Bundle preprocessing for numerical and categorical data
    preprocessor = ColumnTransformer(transformers=[
        ('num', numerical_transformer, numerical_cols),
        ('cat', categorical_transformer, categorical_cols)
    ])

    # Define the model
    model = RandomForestRegressor(n_estimators=100, criterion='squared_error', max_depth=None,
                                  min_samples_split=2, min_samples_leaf=1, max_features=1.0,
                                  bootstrap=True)

    # Bundle preprocessing and modeling code in a pipeline
    my_pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                                  ('model', model)])

    # # Splitting data into training and validation sets
    # X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=0)

    # Preprocessing of training data, fit model
    my_pipeline.fit(X, y)

    return my_pipeline

def predict(trained_model, X_new):
    """
    Makes predictions with the trained model pipeline.

    Parameters:
    - trained_model: The trained model pipeline.
    - X_new: New data to make predictions on.

    Returns:
    - The predictions.
    """
    y_pred = trained_model.predict(X_new)
    return pd.DataFrame(y_pred, columns=['Predictions'])

# Example usage
# Assume df is your DataFrame and 'target' is your target column
# X = df.drop('target', axis=1)
# y = df['target']
# trained_model = train_model(X, y)
# predictions = predict_model(trained_model, X_new) # X_new is your new data for prediction
