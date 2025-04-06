## Project Overview

This repository contains code and data for the paper (Online now) **"TBESO-BP: an improved regression model for predicting subclinical mastitis"**. In this study, we developed a neural network model based on a backpropagation (BP) neural network enhanced with the Multi-strategy Boosted Snake Optimizer (TBESO) to predict somatic cell count (SCC) and assist in the early detection of subclinical mastitis in dairy cows.

The data processing and model evaluation scripts are implemented in Python Jupyter Notebooks, while the machine learning models, including TBESO-BP, alternative BP models, and TBESO ablation experiments, are implemented in MATLAB.

## Repository Contents

- `data_processing/`: Jupyter Notebooks for data processing, cleaning and preparation.
- `models/`: MATLAB scripts for implementing the TBESO-BP model, alternative BP models, and TBESO ablation study experiments, and the code for evaluating models performance metrics calculation such as R², MAE, and RMSE is also included.
- `results/`: Sample results, including model performance metrics and visualizations.
- `evaluation/`: Files and directories for running the TBESO-BP model's optimization and evaluation functions.
- `data/`: Processed dataset used in the study (raw data is not included due to privacy concerns).
- `Population Distribution/`: Visualize the population distribution.

## Getting Started

### Prerequisites

To run the data processing scripts in Python, you will need the following packages:

- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `rcParams`

To run the model training and ablation experiments, you will need **MATLAB** with support for the Neural Network Toolbox (or equivalent) for BP models.

### Directory Structure

- data_processing/

  : Includes:
   - `data_process_of_DHI.ipynb`: Jupyter Notebooks script for data processing, cleaning and preparation.

- models/

  : Includes:

  - `classic_BP.m`: MATLAB script for the classic BP model.
  - `AHL_BP.m`: MATLAB script for the adaptive BP model.
  - `PSO_BP.m`: MATLAB script for the BP model optimized with Particle Swarm Optimization (PSO).
  - `SO_BP.m`: MATLAB script for the BP model optimized with Snake Optimization.
  - `BEESO_BP.m`: MATLAB script for the BP model optimized with BEESO.
  - `TBESO_BP.m`: MATLAB script for the TBESO-enhanced BP neural network model.
  
- evaluation/

  : Includes:

  - `main.m` and `main_data.m`: Main MATLAB scripts to initialize the TBESO-BP model evaluation.
  - `curve_display.m`: Displays optimization curves for evaluation.
  - `cec17_func.cpp` and `cec17_func.mexw64`: Source and compiled files for CEC2017 benchmark functions.
  - `Get_Functions_cec2017.m` and `run_CEC2017.m`: Functions to load and execute CEC2017 benchmarks.
  - Subdirectories like `Apperance`, `CEC2017_curve`, `input_data17`, `My Optimization Algorithms`, and `output_data` for storing visualizations, curve data, input/output data, and additional optimization algorithms.

- results/

  : Includes:
  
  - `10 runs of the model comparison data & figure.xlsx`: Contains data and visualization comparing the performance of different models across 10 independent runs.
  - `Comparison of actual and predicted values.xlsx`: Provides a side-by-side comparison of actual vs. predicted values for SCC, allowing insight into the model's prediction accuracy.
  - `mean_std.xlsx`: Includes the mean and standard deviation for the metrics across different models, helping to assess the stability and variability in model performance.
  - `wilcoxo_rank.xlsx`: Results from the Wilcoxon rank-sum test, providing statistical comparison between models to validate significant differences in performance.

- Population Distribution/
  : includes:

  - `draw.py`: Provides the python to visualize the population distribution.
  - `data/`: each iterations data samples.

### Running the Code

1. **Data Processing**
   Run the data processing scripts in Python to clean and prepare the dataset:

2. **Training the Models**
   To train the models, open each `.m` file in MATLAB (e.g., `tbseobp_model.m` for the TBESO-BP model) and run it. Ensure that your data is loaded correctly within MATLAB as described in the script comments.

3. **Running the TBESO Ablation Study**
   Open `evaluation/main.m` in MATLAB to run ablation tests that evaluate the impact of individual strategies in the TBESO algorithm.

### Data Availability

The dataset used in this project is based on Dairy Herd Improvement (DHI) data. Processed data is available in the `data/` folder. Note that raw data is not included due to privacy restrictions. For more information on the dataset, please refer to the data processing scripts in `data_processing/`.

The dataset contains various production indicators for dairy cows, with detailed descriptions of each variable as follows:

| Variable Name               | Unit    | Description                                                  |
| --------------------------- | ------- | ------------------------------------------------------------ |
| Parity                      | Count   | Number of times the cow has given birth.                     |
| Lactation Days              | Days    | Number of days in the lactation period, starting from the calving date. |
| Milk Yield                  | Kg      | Amount of milk produced by the cow during the lactation period. |
| Milk Fat Rate               | %       | Percentage of fat in the milk.                               |
| Protein Rate                | %       | Percentage of protein in the milk.                           |
| Urea Nitrogen               | mg/dl   | Concentration of urea nitrogen in milk, reflecting the cow’s protein metabolism status. |
| Persistence                 | None    | Lactation persistence, assessing the cow's sustained production ability. |
| Previous Somatic Cell Count | 10^4/ml | Somatic cell count from the previous measurement.            |
| Previous Somatic Cell Score | None    | Somatic cell score from the previous measurement.            |
| Peak Day                    | Days    | The day on which the lactation yield peaked.                 |
| Month                       | Month   | Month in which the data was collected.                       |
| Lactose                     | %       | Percentage of lactose in the milk.                           |
| Somatic Cell Count          | 10^4/ml | Number of somatic cells per milliliter of milk, used for health monitoring. |

## Model Description

The TBESO-BP model leverages a BP neural network framework, enhanced with the TBESO algorithm to optimize initial weights and thresholds. This approach addresses issues of local minima and low population diversity, which are common in traditional BP networks. In addition to the TBESO-BP model, we implemented five alternative BP-based models in MATLAB:

- **Classic BP**: A standard backpropagation neural network.
- **Adaptive BP**: A BP network with adaptive learning rate adjustments.
- **PSO-BP**: A BP network optimized using Particle Swarm Optimization (PSO).
- **Snake-BP**: A BP network optimized using the Snake Optimization algorithm.
- **BEESO-BP**: A BP network optimized with the Boosted Enriched Elite Search Optimization (BEESO) algorithm.

Additionally, the **TBESO Ablation Study** evaluates the individual components of the TBESO algorithm, analyzing their individual contributions to model accuracy.

For detailed algorithmic explanations, please refer to the paper.

## Results

The TBESO-BP model achieved high accuracy in predicting SCC, as demonstrated by the following metrics on the test dataset:

- **R²**: 0.94
- **MAE**: 2.07
- **RMSE**: 5.04

## Contact

For questions or additional information, please contact me at averil_tung@163.com.
