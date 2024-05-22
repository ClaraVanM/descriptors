import pandas as pd
from sklearn.svm import SVC
from skopt.space import Real, Categorical, Integer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from skopt import BayesSearchCV


########################################################################################################################
df1 = pd.read_csv("descriptors_not3.2.1.csv", index_col=0)
df1 = df1[df1["AAC_C"] != -2]
df1["group"] = 0

df2 = pd.read_csv("descriptors_ec3.2.1.csv", index_col=0)
df2 = df2[df2["AAC_C"] != -2]
df2['group'] = 1
df = pd.concat([df1, df2])

nan = df.isna().sum()
nan_df = df[df.isna().any(axis=1)]
nan_columns = df.columns[df.isna().any()]
nan_df = pd.concat([nan_df['name'], nan_df[nan_columns]], axis=1)
nan_df.to_csv("nan.csv")
print(list(nan[nan>0]))



########################################################################################################################

"""X = df.drop(columns=["group", 'name'])
X.to_csv('X.csv')
##sometimes narrowness and corresponding buriedness composition is NaN?
#for now put NaN to 0
X = X.fillna(0)
y = df["group"]
scaler = StandardScaler()
X = scaler.fit_transform(X)

#stratified split because response is unbalances, now the training set has overall an equal amount of both groups as the original dataset
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.75, random_state=0, stratify=y)

#define parameter space
opt = BayesSearchCV(SVC(class_weight='balanced'), {
    'C':Real(0.1, 100, "log-uniform"),
    'gamma':Real(0.01, 10, "log-uniform"),
    'degree': Integer(1,8),
    'kernel': Categorical(["rbf", "linear", "poly"])},
           n_iter=32,
            random_state=0)

_ = opt.fit(X_train, y_train)
print(opt.score(X_test, y_test))"""

### accuracy = 0.90376