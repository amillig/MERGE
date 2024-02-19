# version         v0.1.7
# date            04.02.2024
# author          Alexander-Maurice Illig
# affilation      Chair of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import merge
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from scipy.stats import spearmanr

## User input
startingPosition = 170
paramsFile = 'yap1.params'
csvFile = 'yap1.csv' # ensure that separator is ';' and not ','
fitnessColumn = 'linear'
yWt = 1.0
maxSubstitutions = 3
factor = 1.1
nWalkers = 96
nCores = 8

## Generate sequence representations
print("Initializing 'Encode' class ...")
encodeCls = merge.Encode(startingPosition, paramsFile)
np.save('yap1_wt_encoded.npy', encodeCls._encode_wt())

print('Generating sequence representations ...')
data = merge.get_data(fitnessColumn, csvFile, encodeCls)
merge.generate_dataframe(data, 'yap1_encoded.csv')

## Train model
df = pd.read_csv('yap1_encoded.csv', sep=';')
xWt = np.load('yap1_wt_encoded.npy')
x = df.iloc[:,2:].to_numpy()
deltaE = merge.X_to_deltaE(x, xWt)
y = df.iloc[:,1].to_numpy()   

print('Training model using 80 % of the data ...')
xTrain, xTest, deltaETrain, deltaETest, yTrain, yTest = train_test_split(x, deltaE, y, train_size=0.8, random_state=42) 
model = merge.CombinedPredictor()
model.train(xTrain, deltaETrain, yTrain)

## Predict fitness from sequence
print('R2 score (test set): %.2f '%(r2_score(yTest, model.predict(xTest, deltaETest))))
print('Spearman rho (test set): %.2f '%(spearmanr(yTest, model.predict(xTest, deltaETest))[0]))

## Explore fitness landscape
print('Exploring landscape ...')

explore = merge.Explore(encodeCls, model, yWt, maxSubstitutions=maxSubstitutions, factor=factor)
results = explore.scrape_landscape(nWalkers, nCores=nCores)
df = pd.DataFrame(results, columns=['Variant', 'Predicted Fitness'])
df.to_csv('improved_variants.csv', index=False, sep=';')