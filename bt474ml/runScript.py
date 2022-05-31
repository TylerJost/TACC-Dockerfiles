# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import cross_val_score


from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
# %% Load data
adata = sc.read_h5ad('../data/anndataObjects/BT474LineageAssigned.h5ad')
adata.obs['sample'] = adata.obs['sample'].str.replace('\d+','')
adata = adata[adata.obs['collapsedLineages'] != 'nan']


linFreqs = adata.obs.groupby(['collapsedLineages', 'sample']).size().unstack().reset_index().rename_axis(None, axis=1).set_index('collapsedLineages')
# Add 1 count to everything and convert to proportions
linFreqs = linFreqs+1
linFreqs = linFreqs/np.sum(linFreqs, axis=0)
linFreqs.head()
# %% Convert to fold changes
def logFC(F,O):
    """
    Returns the log2 fold change
    Take final proportion and initial proportion
    """
    return np.log2(F/O)
comparisons = [['D', 'PreTreat'], ['LP','PreTreat'], ['DLP','D'], ['LPD','LP']]

fcs = []
fcNames = []
for comparison in comparisons:
    finalSample = linFreqs[comparison[0]]
    initialSample = linFreqs[comparison[1]]
    fcs.append(logFC(finalSample, initialSample))
    fcNames.append(comparison[0]+'-FC')
fcs = pd.DataFrame(fcs).transpose()
fcs.columns = fcNames
fcs = fcs.reset_index()
fcs.head()

adata.obs = adata.obs.merge(fcs, how='left').set_index(adata.obs.index)
adataPre = adata[adata.obs['sample'] == 'PreTreat',]
# %% SVM Regression
print('Running SVM!')
X = adataPre.X
varRegress = 'D-FC'
y = list(adataPre.obs[varRegress])

svr = GridSearchCV(
    SVR(kernel="rbf", gamma=0.1),
    param_grid={"C": [1e0, 1e1, 1e2, 1e3], "gamma": np.logspace(-2, 2, 5)},
)
svr.fit(X, y)

print(f"Best SVR with params: {svr.best_params_} and R2 score: {svr.best_score_:.3f}")
dfScores = pd.DataFrame(svr.cv_results_)
dfScores.to_csv('../data/scoresSVM')
# %%

