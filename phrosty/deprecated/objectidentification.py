# ruff: noqa
# (Skip ruff, this is deprecated anyway!)

# IMPORTS Standard:
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, confusion_matrix, ConfusionMatrixDisplay, recall_score
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.model_selection import train_test_split
import pickle as pkl
import pandas as pd
import numpy as np
import os
from operator import methodcaller

# IMPORTS Astro: 
from astropy.table import Table

# IMPORTS Internal: 
from .utils import train_config, predict_config
from .plotting import classification_contours

"""
Train a Random Forest Classifier on star/galaxy truth coordinates.
This model is then used in photometry.py to generate the ePSF. 
"""

def train_model(train_table,test_size=0.4,random_state=42):
    """
    Inputs: flux, max pixel, ellipticity for galaxies and stars. 
    Astropy table with sum of flux in pixels, max pixel value, ellipticity, and object type.

    """

    train_table = train_table.to_pandas()
    X = train_table[['peak','flux','ellipticity']]
    y = np.array(list(map(train_config,train_table['objtype'])))

    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=test_size,random_state=random_state)

    model = RandomForestClassifier(max_depth=100, n_estimators=100, max_features=1)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    score = model.score(X_test, y_test)
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average='weighted')
    recall = recall_score(y_test, y_pred, average='weighted')

    data = {'X': X, 'y': y, 
            'X_train': X_train, 'X_test': X_test, 
            'y_train': y_train, 'y_test': y_test,
            'y_pred': y_pred}
    scores = {'score': score, 'accuracy': accuracy, 'precision': precision, 'recall': recall}

    return model, data, scores

def plot_contours(data, model, **kwargs):
    classification_contours(data, model, **kwargs)

def plot_confusionmatrix(data):
    cm = confusion_matrix(data['y_test'], data['y_pred'])
    ConfusionMatrixDisplay(confusion_matrix=cm,display_labels=['Galaxy', 'Star']).plot()

def classify(predict_table, modelname):
    roman_bands = ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'W146', 'K213']
    if modelname in roman_bands:
        cdir = os.path.dirname(os.path.abspath(__file__))
        modelpath = os.path.join(cdir,f'models/{modelname}_objid_model.pkl')
    else:
        modelpath = copy(modelname)
    
    modelfile = open(modelpath, 'rb')
    model = pkl.load(modelfile)

    predict_table = predict_table.to_pandas()
    X = predict_table[['peak','flux','ellipticity']]
    X_shape = X.shape
    boolean_array = np.empty(X_shape)

    for i, col in enumerate(X):
        boolean_array[:,i] = np.isnan(X[col])

    method_any = methodcaller('any')
    nanshere = list(map(method_any, boolean_array))
    nonanshere = ~np.array(nanshere)
    nonansidx = np.where(nonanshere)[0]

    predicted = model.predict(X.iloc[nonansidx])
    pred_type = np.array(list(map(predict_config,predicted)))

    predict_table = Table.from_pandas(predict_table)
    predict_table['predicted objtype'] = 'other'
    predict_table['predicted objtype'][nonansidx] = pred_type
    
    return predict_table
