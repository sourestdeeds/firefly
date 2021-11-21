from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, MinMaxScaler, StandardScaler, RobustScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from xgboost.sklearn import XGBRegressor
from xgboost import plot_tree, plot_importance, to_graphviz
from sklearn.model_selection import GridSearchCV, KFold
from sklearn import preprocessing
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from typing import List
from pathlib import Path



def read_data(filename: str, objective: str, features: List[str], seed=10, override_features=False):
    # Read the data
    global X_train, X_test, y_train, y_test, X, y
    # Add the objective to cols
    featuresAndobjective = features.copy()
    featuresAndobjective += [objective]
    filename = Path(filename)
    if override_features:
        X = pd.read_csv(filename)
    else:
        X = pd.read_csv(filename, usecols=featuresAndobjective)

    # Remove rows with missing target, separate target from predictors
    X.dropna(axis=0, subset=[objective, 'pl_radj'], inplace=True)
    # X.dropna(axis=0, subset=[objective], inplace=True)
    y = X[objective]
    # print(f'Target has {y.isnull().sum()} missing values.')
    X.drop([objective], axis=1, inplace=True)

    # Break off validation set from training data
    X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                            train_size=0.8, test_size=0.2,
                                            random_state=seed)
    #print(X_train.columns)
    return X_train, X_test, y_train, y_test, X, y


def pipelines(hyper_parameters, seed=0):
    '''
    Runs the preprocessing and model pipeline.
    '''

    # "Cardinality" means the number of unique values in a column
    # Select categorical columns with relatively low cardinality (convenient but arbitrary)

    categorical_cols = [cname for cname in X_train.select_dtypes('object').columns 
                        if X_train[cname].nunique() < 10] 

    # Select numerical columns
    numerical_cols = X_train.select_dtypes('number').columns 

    # Preprocessing for numerical data
    numerical_transformer =  Pipeline(
        steps=[
            ('imputer', SimpleImputer(strategy='mean')),
            ('minmaxscaler', RobustScaler()),
        ])

    # Preprocessing for categorical data
    categorical_transformer = Pipeline(
        steps=[
            ('imputer', SimpleImputer(strategy='most_frequent')),
            ('onehot', OneHotEncoder(handle_unknown='ignore', sparse=False)),
        ])

    # Bundle preprocessing for numerical and categorical data
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', numerical_transformer, numerical_cols),
            ('cat', categorical_transformer, categorical_cols)
        ])

    cv = KFold(n_splits=5) 
    # Tune Params and cross validate
    xgbRegressor_gridSearch = Pipeline(
        steps=[
            ('crossValidation', GridSearchCV(
                XGBRegressor(random_state=seed),
                hyper_parameters,
                cv = cv,
                n_jobs = 10,
                verbose=True,
                refit=True
                )
            )
        ])
    
    # Run the entire pipeline
    clf = Pipeline(
        steps=[
            ('preprocessor', preprocessor),
            ('model', xgbRegressor_gridSearch)
        ])
    
    # Preprocessing of training data, fit model 
    clf.fit(X_train, y_train)
    # Refit with the best params
    # Preprocessing of validation data, get predictions
    y_pred = clf.predict(X_test)
    # plot_tree(clf.named_steps['model'])
    # Save the model
    import pickle
    pickle.dump(clf, open('model.pkl', 'wb'))

    print(f'MAE: {mean_absolute_error(y_test, y_pred):.4f}')
    print(f'Accuracy: {clf.score(X_test, y_test) * 100:.2f}%')
    return clf, y_pred

def permutations():
    import eli5
    from eli5.sklearn import PermutationImportance


    # Break off validation set from training data numeric only
    X_train, X_test, y_train, y_test = train_test_split(X.select_dtypes('number'), y, 
                                                                    train_size=0.8, test_size=0.2,
                                                                    random_state=1)

    categorical_cols = [cname for cname in X_train.select_dtypes('object').columns 
                            if X_train[cname].nunique() < 10] 

    # Select numerical columns
    numerical_cols = X_train.select_dtypes('number').columns 

    # Preprocessing for numerical data
    numerical_transformer =  Pipeline(
        steps=[
            ('imputer', SimpleImputer(strategy='mean')),
            ('minmaxscaler', MinMaxScaler(feature_range=(0,1)))
        ])

    # Preprocessing for categorical data
    categorical_transformer = Pipeline(
        steps=[
            ('imputer', SimpleImputer(strategy='most_frequent')),
            ('onehot', OneHotEncoder(handle_unknown='ignore', sparse=False)),
        ])

    # Bundle preprocessing for numerical and categorical data
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', numerical_transformer, numerical_cols),
            ('cat', categorical_transformer, categorical_cols)
        ])

    # Run the entire pipeline
    permutation_test = Pipeline(
        steps=[
            ('preprocessor', preprocessor),
            ('model', XGBRegressor(n_estimators=500, random_state=0))
        ])


    permutation_test.fit(X_train, y_train)

    perm = PermutationImportance(permutation_test, random_state=0).fit(X_test, y_test)
    #eli5.show_weights(perm, feature_names = X_test.drop('pl_name', axis=1).columns.tolist())
    return perm


def compare_predict_vs_test(raw_data, X_test, y_test, y_pred, mass=True, sample_size=15):
    '''Prints a table to compare the test versus the prediction.'''
    from tabulate import tabulate
    import mr_forecast as mr
    
    
    # Prediction
    y_pred = pd.Series(y_pred).rename('Prediction')
    y_test = y_test.rename('Test')

    # Set same indexes
    y_pred.index = y_test.index
    
    # Test
    y_test_p = y_test.sample(sample_size)
    X_test_p = X_test['pl_radj'][y_test_p.index]
    y_pred_p = y_pred[y_test_p.index]
    pl_name = raw_data.pl_name[y_test_p.index]

    df = pd.DataFrame()
    df['Planet'] = pl_name
    df['Test'] = y_test_p
    df['Prediction'] = y_pred_p
    if mass:
        print('Predicting the Mass:\n')
        df['Forecaster'] = [mr.Rstat2M(mean=i, std=0.01, unit='Jupiter', sample_size=100, grid_size=100)[0] for i in list(X_test_p)]
    else: # Radius
        print('Predicting the Radius:\n')
        df['Forecaster'] = [mr.Mstat2R(mean=i, std=0.01, unit='Jupiter', sample_size=100, grid_size=100)[0] for i in list(X_test_p)]
    df['Residual Test-Prediction'] = np.abs(df['Test'] - df['Prediction'])
    df['Residual Test-Forecaster'] = np.abs(df['Test'] - df['Forecaster'])
    df['Prediction < Forecaster'] = df['Residual Test-Prediction'] < df['Residual Test-Forecaster']
    print(tabulate(df, tablefmt='psql', showindex=False, headers='keys', floatfmt=".5f"))
    print(f"\nResidual sum for Prediction: {df['Residual Test-Prediction'].sum():.2f}\nResidual sum for Forecaster: {df['Residual Test-Forecaster'].sum():.2f}")
    print(f"Prediction versus Forecaster Accuracy: {df['Prediction < Forecaster'].sum()}/{sample_size}")
    return df



if __name__ == '__main__':
    # Parameter Tuning
    hyper_parameters = {'nthread':[4], #when use hyperthread, xgboost may become slower
                        'objective':'count:poisson',#['reg:squarederror'],
                        'learning_rate': [.04], #so called `eta` value
                        'max_depth': [3],
                        'min_child_weight': [0],
                        'subsample': [0.2],
                        'colsample_bytree': [1],
                        'n_estimators': [500]}


    X_train, X_test, y_train, y_test, X, y = read_data('data/ML_nasa_tess_viable.csv.xz', 'pl_massj', seed=3)
    #X_train, X_test, y_train, y_test, X, y = read_data('data/nasa_full.csv.xz', 'pl_massj', seed=3)
    #X_train, X_test, y_train, y_test, X, y = read_data('data/nasa_product_list.csv.xz', 'pl_bmassj')
    #X_train, X_test, y_train, y_test, X, y = read_data('data/kep.csv.gz', 'pl_massj')
    # usecols=['pl_radj', 'pl_massj']
    clf, y_pred = pipelines(hyper_parameters, seed=0)


    mass_cols = ['pl_cmassj',
            'pl_cmassjerr1',
            'pl_cmassjerr2',
            'pl_cmassjlim',
            'pl_cmassjstr',
            'pl_cmasse',
            'pl_cmasseerr1',
            'pl_cmasseerr2',
            'pl_cmasselim',
            'pl_cmassestr',
            #'pl_massj',
            'pl_massjerr1',
            'pl_massjerr2',
            'pl_massjlim',
            'pl_massjstr',
            'pl_masse',
            'pl_masseerr1',
            'pl_masseerr2',
            'pl_masselim',
            'pl_massestr',
            'pl_bmassjerr1',
            'pl_bmassjerr2',
            'pl_bmassjlim',
            'pl_bmassjstr',
            'pl_bmasseerr1',
            'pl_bmasseerr2',
            'pl_bmasselim',
            'pl_bmassestr',
            'pl_bmassprov',
            'pl_msinij',
            'pl_msinijerr1',
            'pl_msinijerr2',
            'pl_msinijlim',
            'pl_msinijstr',
            'pl_msinie',
            'pl_msinieerr1',
            'pl_msinieerr2',
            'pl_msinielim',
            'pl_msiniestr',
            'pl_bmassj',
            'pl_bmasse',
            'pl_dens', 
            'pl_denserr1',
            'pl_denserr2', 
            'pl_denslim', 
            'pl_densstr',
            'pl_rvamp',
            'pl_rvamp', 
            'pl_rvamperr1', 
            'pl_rvamperr2', 
            'pl_rvamplim',
            'pl_rvampstr']

    radius_cols = [#'pl_radj',
                'pl_radjerr1',
                'pl_radjerr2',
                'pl_radjlim',
                'pl_radjstr',
                'pl_rade',
                'pl_radeerr1',
                'pl_radeerr2',
                'pl_radelim',
                'pl_radestr',]