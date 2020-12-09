""" Post processes individual experiment results into a single csv file for analysis and visualization """

import pandas as pd
import os

experiment_results_file = 'experiment_results/'
df = pd.concat([pd.read_csv(experiment_results_file + file) for file in os.listdir(experiment_results_file)])
df.index = os.listdir(experiment_results_file)
df.to_csv('experiment_results.csv')
