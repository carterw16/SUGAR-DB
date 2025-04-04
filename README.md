# Microgrid Simulation Tool

An interactive web-based simulation tool for designing and analyzing microgrid systems. This tool combines real-time visualization, machine learning-based forecasting, and power flow optimization to support research, education, and practical use in microgrid planning and operations.

## Features

- **Custom Microgrid Modeling** – Upload your own GridLAB-D-based microgrid designs and simulate node behavior.
- **ML Forecasting** – Predict solar, wind, and load data using RNN-based models trained on historical weather datasets.
- **Power Flow Visualization** – Real-time interactive visualization of power flow, node status, and load balancing using Vis.js.
- **Weather-Driven Forecasts** – Integrates weather data from NSRDB and Kaggle to drive energy prediction models.
- **Optimization** – Multi-period AC-OPF with DDP-based convergence using the SUGAR3 framework.
- **Web Application** – Built with Django, Bootstrap, and deployed on AWS EC2 (currently disabled) for accessible and scalable performance.

## Tech Stack

- **Backend:** Django, Python, SUGAR3
- **Frontend:** HTML/CSS, Bootstrap, Vis.js
- **Machine Learning:** LSTM/RNN, scikit-learn, NumPy, Pandas
- **Deployment:** AWS EC2, SQLite
- **Data Sources:** NSRDB, Kaggle Wind Turbine Data

## Getting Started

To run the project locally:

1. Clone this repo:
   ```bash
     git clone https://github.com/your-username/microgrid-simulator.git
     cd microgrid-simulator
   ```
2. Create and activate virtual environment and install dependencies:
   ```bash
    python -m venv venv
    source venv/bin/activate  # or `venv\Scripts\activate` on Windows
    pip install -r requirements.txt
   ```
3. Run the server
   ```bash
    pip install -r requirements.txt
   ```
4. Visit `http://localhost:8000`in your browser to access the app.
