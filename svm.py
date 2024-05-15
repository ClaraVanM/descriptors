import pandas as pd
from sklearn.svm import SVC
from skopt.space import Real, Categorical, Integer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from skopt import BayesSearchCV


########################################################################################################################
df1 = pd.read_csv("descriptors_not3.2.1_new.csv", index_col=0)
df1 = df1[df1["AAC_C"] != "-2"]
df1["group"] = 0

df2 = pd.read_csv("descriptors_ec3.2.1_new.csv", index_col=0)
df2 = df2[df2["AAC_C"] != "-2"]
df2['group'] = 1
df = pd.concat([df1, df2])
structure_based = ['depth', 'narrowness_0', 'narrowness_1', 'narrowness_2', 'narrowness_3','buried_CYS1','buried_ASP1','buried_SER1','buried_GLN1','buried_LYS1','buried_ILE1','buried_PRO1','buried_THR1','buried_PHE1','buried_ASN1','buried_GLY1','buried_HIS1','buried_LEU1','buried_ARG1','buried_TRP1','buried_ALA1','buried_VAL1','buried_GLU1','buried_TYR1','buried_MET1','buried_CYS2','buried_ASP2','buried_SER2','buried_GLN2','buried_LYS2','buried_ILE2','buried_PRO2','buried_THR2','buried_PHE2','buried_ASN2','buried_GLY2','buried_HIS2','buried_LEU2','buried_ARG2','buried_TRP2','buried_ALA2','buried_VAL2','buried_GLU2','buried_TYR2','buried_MET2','buried_CYS0','buried_ASP0','buried_SER0','buried_GLN0','buried_LYS0','buried_ILE0','buried_PRO0','buried_THR0','buried_PHE0','buried_ASN0','buried_GLY0','buried_HIS0','buried_LEU0','buried_ARG0','buried_TRP0','buried_ALA0','buried_VAL0','buried_GLU0','buried_TYR0','buried_MET0','buried_CYS3','buried_ASP3','buried_SER3','buried_GLN3','buried_LYS3','buried_ILE3','buried_PRO3','buried_THR3','buried_PHE3','buried_ASN3','buried_GLY3','buried_HIS3','buried_LEU3','buried_ARG3','buried_TRP3','buried_ALA3','buried_VAL3','buried_GLU3','buried_TYR3','buried_MET3','buried_CYS4','buried_ASP4','buried_SER4','buried_GLN4','buried_LYS4','buried_ILE4','buried_PRO4','buried_THR4','buried_PHE4','buried_ASN4','buried_GLY4','buried_HIS4','buried_LEU4','buried_ARG4','buried_TRP4','buried_ALA4','buried_VAL4','buried_GLU4','buried_TYR4','buried_MET4','O','C','N','S','H']
df = df.drop(columns=structure_based)
"""nan = df.isna().sum()
print(nan[nan>0])"""



########################################################################################################################

X = df.drop(columns=["group", 'name'])
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
opt = BayesSearchCV(SVC(), {
    'C':Real(0.1, 100, "log-uniform"),
    'gamma':Real(0.01, 10, "log-uniform"),
    'degree': Integer(1,8),
    'kernel': Categorical(["rbf", "linear", "poly"])},
           n_iter=32,
            random_state=0)

_ = opt.fit(X_train, y_train)
print(opt.score(X_test, y_test))

### accuracy = 0.8677