#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from time import time
import re

import numpy as np
import pandas as pd
import scipy
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.grid_search import RandomizedSearchCV, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from scipy.stats import randint as sp_randint
from sklearn.externals import joblib

import matplotlib as mpl 
mpl.use("agg")
from matplotlib import pyplot as plt 
mpl.rc('font', family='serif') 
mpl.rc('text', usetex='true') 
mpl.rc('text', dvipnghack='true') 
mpl.rcParams.update({'font.size': 16}) 
import pylab as P

y_latex = ("$M_0$", "$Y_0$", "$Z_0$", "$\\alpha_{\mathrm{MLT}}$", 
    "\mathrm{Age}")
y_names = ["M", "Y", "Z", "alpha", "age"]

"""
X_names = ["$\mathrm{radius}$", "$L$", "$T_{\mathrm{eff}}$", "$\log g$", 
    "$\mathrm{Fe/H}$", "$\\nu_\max$", 
    "$\Delta \\nu_0$", "$\epsilon_0$", "$\delta \\nu_{02}$", 
    "$\Delta \\nu_1$", "$\epsilon_1$", "$\delta \\nu_{13}$", 
    "$\Delta \\nu_2$", "$\epsilon_2$", 
    "$\Delta \\nu_3$", "$\epsilon_3$"]
"""

### Load data
#data = np.loadtxt('grids/sobol_grid-diffusion.dat', skiprows=1)
data = pd.read_csv('grids/deleter.dat', sep='\t')
exclude_pattern = "He|log_g|mass|nu_max|Dnu"
data = data.drop([i for i in data.columns 
    if re.search(exclude_pattern, i)], axis=1)
X = data.drop(y_names, axis=1)
ys = data.drop(X.columns, axis=1)
#ys['age'] = np.log10(ys['age']/10**9)# / 10**9
print(X.head())
#X = np.array(data)[:,5:]
#del_indices = [X_names.index("Dnu_3"), X_names.index("eps_3")]
#X = np.delete(X, del_indices, 1)
#X_names = np.delete(X_names, del_indices)
#ys = np.array(data)[:,:5]
#ys[:,-2] = ys[:,-2]/10**9

### Train model
"""
param_dist = {"max_features": sp_randint(1, len(X_names)),
              "min_samples_split": sp_randint(1, len(X_names)),
              "min_samples_leaf": sp_randint(1, len(X_names))}
forest = RandomizedSearchCV(
    ExtraTreesRegressor(1000, n_jobs=-1, verbose=1), 
    param_distributions=param_dist, n_iter=100, n_jobs=-1, verbose=1)
"""
"""
forest = GridSearchCV(
    RandomForestRegressor(1000, n_jobs=-1, verbose=1, bootstrap=1, oob_score=1),
    param_grid={"min_samples_split": list(range(1, len(X_names)))},
    verbose=1, n_jobs=-1)
"""
"""
y_trfm = Pipeline(steps=[
    ('scaler', StandardScaler()), 
    ('pca', PCA())
])
new_ys = y_trfm.fit_transform(ys)

forest = Pipeline(steps=[
    ('scaler', StandardScaler()), 
    ('pca', PCA()),
    ('rf', RandomForestRegressor(n_estimators=100, n_jobs=1, verbose=1,
            oob_score=1, bootstrap=1))])
"""
forest = ExtraTreesRegressor(n_estimators=2000, n_jobs=-1, verbose=1,
    oob_score=1, bootstrap=1)
start = time()
forest.fit(X, ys)#new_ys)
end = time()
print("%.5g seconds to train regressor" % (end-start))
est = forest#.steps[2][1]#.best_estimator_
print(est.oob_score_)
#joblib.dump(forest, 'rf.pickle')

### Run on all stars in the perturb directory
print("Star\t\tM\t\tY\t\tZ\t\talpha\t\tage")
for star_fname in os.listdir('perturb'):
    #star_data = np.loadtxt(os.path.join('perturb', star_fname), skiprows=1)
    #star_data = np.delete(star_data, del_indices, 1)
    star_data = pd.read_csv(os.path.join('perturb', star_fname), sep='\t')
    star_data = star_data.drop([i for i in star_data.columns 
        if re.search(exclude_pattern, i)], axis=1)
    predict = forest.predict(star_data)
    #predict = y_trfm.inverse_transform(predict)
    #predict[:,y_names.index('age')] /= 10**9
    #predict[:,y_names.index('age')] = 10**predict[:,y_names.index('age')]
    middles = np.median(predict, 0)
    stds = np.std(predict, 0)
    outstr = star_fname.split("_")[0]
    plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.suptitle(outstr)
    for (col_j, name) in enumerate(y_names):
        (m, s) = (middles[col_j], stds[col_j])
        outstr += "\t%.3g +/- %.3g" % (m, s)
        if (col_j == 4):
            ax = plt.subplot2grid((3,4), (2,1), colspan=2)
        else:
            ax = plt.subplot2grid((3,4), (int(col_j/2),(col_j%2)*2), colspan=2)
        #ax = plt.subplot(3,2,col_j+1)
        
        n, bins, patches = P.hist(predict[:,col_j], 50, normed=1, 
            histtype='stepfilled', color='white')
        #y = P.normpdf(bins, m, s)
        #P.plot(bins, y, 'b--', linewidth=1.5)
        
        P.xlabel(y_latex[col_j])
        ax.set_yticklabels('',visible=False)
        P.locator_params(axis='x', nbins=5)
        xticks = [m-4*s, m-2*s, m, m+2*s, m+4*s]
        ax.set_xlim([m-4*s, m+4*s])
        ax.set_xticks(xticks)
        ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
        ax.minorticks_on()
        plt.tight_layout()
    
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join('learn_plots', star_fname.split("_")[0] + '.png'))
    """
    plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.matshow(pd.DataFrame(predict, columns=y_names).corr())
    plt.savefig(os.path.join('learn_plots', 
        star_fname.split("_")[0] + '-corr.png'))
    """
    print(outstr)

### Test on 16 Cyg A & B
#cygA = np.loadtxt('perturb/16CygA_perturb.dat', skiprows=1)
#cygB = np.loadtxt('perturb/16CygB_perturb.dat', skiprows=1)
#cygA = np.delete(cygA, del_indices, 1)
#cygB = np.delete(cygB, del_indices, 1)
cygA = pd.read_csv('perturb/16CygA_perturb.dat', sep='\t')
cygA = cygA.drop([i for i in cygA.columns 
    if re.search(exclude_pattern, i)], axis=1)
cygB = pd.read_csv('perturb/16CygB_perturb.dat', sep='\t')
cygB = cygB.drop([i for i in cygB.columns 
    if re.search(exclude_pattern, i)], axis=1)
start = time()
#pred_A = y_trfm.inverse_transform(forest.predict(cygA))
#pred_B = y_trfm.inverse_transform(forest.predict(cygB))
pred_A = forest.predict(cygA)
pred_B = forest.predict(cygB)
"""
pred_A = y_trfm.inverse_transform(pred_A)
pred_A[:,y_names.index('age')] /= 10**9
pred_B = y_trfm.inverse_transform(pred_B)
pred_B[:,y_names.index('age')] /= 10**9
"""
end = time()
print("%.3g seconds to predict Cyg A and B" % (end-start))
plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
plt.suptitle("16 Cyg")
for (col_j, name) in enumerate(y_names):
    if (col_j == 4):
        ax = plt.subplot2grid((3,4), (2,1), colspan=2)
    else:
        ax = plt.subplot2grid((3,4), (int(col_j/2),(col_j%2)*2), colspan=2)
    #ax = plt.subplot(3,2,col_j+1)
    
    (m_A, s_A) = (np.median(pred_A[:,col_j]), np.std(pred_A[:,col_j]))
    (m_B, s_B) = (np.median(pred_B[:,col_j]), np.std(pred_B[:,col_j]))
    
    n, bins, patches = P.hist(pred_A[:,col_j], 50, normed=1, 
        histtype='stepfilled', color='red', alpha=0.5)
    #y = P.normpdf(bins, m_A, s_A)
    #P.plot(bins, y, 'r--', linewidth=1.5)
    
    n, bins, patches = P.hist(pred_B[:,col_j], 50, normed=1, 
        histtype='stepfilled', color='blue', alpha=0.5)
    #y = P.normpdf(bins, m_B, s_B)
    #P.plot(bins, y, 'b--', linewidth=1.5)
    
    P.xlabel(y_latex[col_j])
    ax.set_yticklabels('',visible=False)
    P.locator_params(axis='x', nbins=5)
    xticks = np.linspace(min(m_A-4*s_A, m_B-4*s_B), 
                         max(m_A+4*s_A, m_B+4*s_B), 5)
    ax.set_xlim([min(xticks), max(xticks)])
    ax.set_xticks(xticks)
    ax.set_xticklabels(['%.3g'%xtick for xtick in xticks])
    ax.minorticks_on()
    plt.tight_layout()

plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join('learn_plots', '16Cyg.png'))

"""
print("\t16CygA\t16CygB")
for (col_j, name) in enumerate(y_names):
    cyg_A_mean = np.mean(pred_A[:,col_j])
    cyg_A_std = np.std(pred_A[:,col_j])
    cyg_B_mean = np.mean(pred_B[:,col_j])
    cyg_B_std = np.std(pred_B[:,col_j])
    if name=="age":
        cyg_A_mean = cyg_A_mean / 10**9
        cyg_A_std  = cyg_A_std  / 10**9
        cyg_B_mean = cyg_B_mean / 10**9
        cyg_B_std  = cyg_B_std  / 10**9
    print("%s\t%.3g ± %.3g\t%.3g ± %.3g" % (name, 
            cyg_A_mean, cyg_A_std, 
            cyg_B_mean, cyg_B_std))
"""

### Plot importances 
#est = forest.steps[1][1]
importances = est.feature_importances_
indices = np.argsort(importances)[::-1]
import_dist = np.array([tree.feature_importances_ 
    for tree in est.estimators_]).T[indices][::-1].T

X_names = X.columns
print("Feature ranking:")
for f in range(len(X_names)):
    print("%d. %s (%.3g)" % (f+1, X_names[indices[f]], 
        importances[indices[f]]))

mpl.rc('text', usetex='false') 
plt.figure(figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
plt.boxplot(import_dist, vert=0)
plt.yticks(range(1,1+len(X_names)), np.array(X_names)[indices][::-1])
plt.xlabel("Feature importance")
plt.tight_layout()
plt.savefig(os.path.join('learn_plots', 'feature_importance.pdf'))
plt.close()
