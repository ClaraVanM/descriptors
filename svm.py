import pandas as pd
from sklearn import  svm

"""df1 = pd.read_csv("descriptors_not3.2.1.csv", index_col=0)
df1 = df1[df1["name"] != "-2"]"""

from sklearn.datasets import load_breast_cancer
from sklearn.svm import SVC
from skopt import gp_minimize
from skopt.space import Real, Categorical
from skopt.utils import use_named_args
from sklearn.model_selection import cross_val_score

data = load_breast_cancer()
X,y = data.data, data.target

#define parameter space
param_space = [Real(0.1, 100, "log-uniform", name = "c"),
               Real(0.01, 10, "log-uniform", name="gamma"),
               Categorical(["rbf", "linear", "poly"], name="kernel")]

#define objective function
@use_named_args(param_space)
def objective(**params):
    clf = SVC(**params)
    accuracy = cross_val_score(clf, X,y, cv=5,n_jobs=1).mean()
    return -accuracy

res_gp = gp_minimize(objective, param_space, n_calls=20, random_state=0)

best_params = dict(zip(['C', 'gamma', 'kernel'], res_gp.x))

print('Best parameters:', best_params)
print('Best accuracy', -res_gp.fun)
