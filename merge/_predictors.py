# version         0.1.7
# date            30.01.2024
# author          Alexander-Maurice Illig
# affilation      Chair of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import numpy as np
from scipy.optimize import curve_fit
class PredictorDCA:
    def __init__(self, maxfev=10000, p0=(1,1,-7,1),
                 bounds=[(-5, -5, -20, -20), (5, 5, 0, 20)]):
        self.maxfev = maxfev
        self.p0 = p0
        self.bounds = bounds
        
    @staticmethod
    def logistic(deltaE, *args):
        return args[0] / (1 + np.exp(-args[1]*(deltaE-args[2]))) + args[3]
    
    def fit(self, deltaE, y):
        self.params, self.cov = curve_fit(self.logistic, deltaE, y, maxfev=self.maxfev, p0=self.p0, bounds=self.bounds)
        return self
    
    def predict(self, deltaE):
        return self.logistic(deltaE, *self.params)
    
    
from sklearn.linear_model import RidgeCV
class PredictorRidge:
    def __init__(self, alphas=np.logspace(-6,6,100),
                 fitIntercept=True, cv=5):
        self.alphas = alphas
        self.fitIntercept = fitIntercept
        self.cv = cv       
    
    def fit(self, x, y):
        self.predictor = RidgeCV(alphas=self.alphas, fit_intercept=self.fitIntercept, cv=self.cv).fit(x, y)
        return self
    
    def predict(self, x):
        return self.predictor.predict(x)


from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
class PredictorLasso:
    def __init__(self, alphas=np.logspace(-6,6,100),
                 fitIntercept=True, cv=5, nJobs=1):
        self.alphas = alphas
        self.fitIntercept = fitIntercept
        self.cv = cv
        self.nJobs = nJobs
    
    def fit(self, x, y):
        grid = GridSearchCV(Lasso(
            fit_intercept=self.fitIntercept),
            {'alpha':self.alphas},
            cv=self.cv,
            n_jobs=self.nJobs
            )
        grid.fit(x, y)
        self.predictor = Lasso(**grid.best_params_, fit_intercept=self.fitIntercept).fit(x, y)
        return self
    
    def predict(self, x):
        return self.predictor.predict(x)


from sklearn.linear_model import LinearRegression
class PredictorOLS:
    def __init__(self, fitIntercept=True):
        self.fitIntercept = fitIntercept
    
    def fit(self, x, y):
        self.predictor = LinearRegression(fit_intercept=self.fitIntercept).fit(x, y)
        return self
    
    def predict(self, x):
        return self.predictor.predict(x)


from sklearn.ensemble import RandomForestRegressor
class PredictorRF:
    def __init__(self,
     nEstimators=[1, 5, 10, 20, 50, 100, 200, 500, 1000],
     maxFeatures=['auto', 'sqrt', 'log2'],
     cv=5, nJobs=1
     ):
        self.nEstimators = nEstimators
        self.maxFeatures = maxFeatures
        self.cv = cv 
        self.nJobs = nJobs

    def fit(self, x, y):
        grid=GridSearchCV(
            RandomForestRegressor(),
            {'n_estimators':self.nEstimators,'max_features':self.maxFeatures},
            cv=self.cv,
            n_jobs=self.nJobs
            )
        grid.fit(x, y)

        self.predictor = RandomForestRegressor(**grid.best_params_).fit(x, y)
        return self
    
    def predict(self, x):
        return self.predictor.predict(x)


from sklearn.svm import SVR
class PredictorSVR:
    def __init__(self,
     Cs=[0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0],
     epsilons=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0],
     cv=5, 
     nJobs=1
     ):
        self.Cs = Cs
        self.epsilons = epsilons
        self.cv = cv 
        self.nJobs = nJobs

    def fit(self, xScaled, y):
        grid=GridSearchCV(
            SVR(),
            {'C':self.Cs, 'epsilon':self.epsilons},
            cv=self.cv,
            n_jobs=self.nJobs
            )
        grid.fit(xScaled, y)

        self.predictor = SVR(**grid.best_params_).fit(xScaled, y)
        return self
    
    def predict(self, xScaled):
        return self.predictor.predict(xScaled)


from sklearn.model_selection import KFold
from scipy.optimize import differential_evolution
class CombinedPredictor:
    def __init__(self, 
        nSplits=5, 
        shuffle=True,
        randomState=12,
        bounds=[(0,10),(0,10)],
        tol=1e-4,
        predictor2=PredictorRidge()
        ):
        self.nSplits = nSplits
        self.shuffle = shuffle
        self.randomState = randomState
        self.bounds = bounds
        self.tol = tol
        self.predictor2 = predictor2

    def five_fold_split(self, x):
        fiveFold = KFold(n_splits=self.nSplits, random_state=self.randomState, shuffle=self.shuffle)
        return fiveFold.split(x)
    
    def train(self, x, deltaE, y):
        data = [[] for _ in range(3)]
        for trainingIdxs,validationIdxs in self.five_fold_split(x):
            data[0].append(y[validationIdxs])
            
            p1 = PredictorDCA()
            p1.fit(deltaE[trainingIdxs],y[trainingIdxs])
            data[1].append(p1.predict(deltaE[validationIdxs]))
            
            p2 = PredictorRidge()
            p2.fit(x[trainingIdxs],y[trainingIdxs])
            data[2].append(p2.predict(x[validationIdxs]))
            
        yTrue,yP1,yP2 = [np.concatenate(l) for l in data]
        
        loss = lambda params: np.sum(np.power(yTrue - params[0]*yP1 - params[1]*yP2, 2))
        minimizer = differential_evolution(loss, bounds=self.bounds, tol=self.tol)        
        self.gamma1, self.gamma2 = minimizer.x
        
        self.p1 = PredictorDCA().fit(deltaE, y)
        self.p2 = self.predictor2.fit(x, y)
        return self
    
    def predict(self, x, deltaE):
        return self.gamma1 * self.p1.predict(deltaE) + self.gamma2 * self.p2.predict(x)
