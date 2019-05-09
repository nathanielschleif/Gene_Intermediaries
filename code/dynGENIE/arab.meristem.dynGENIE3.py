import numpy as np
import pandas as pd
from dynGENIE3 import *
import _pickle
f = open('TS_data.pkl','rb')
(TS_data, time_points, decay_rates, gene_names) = _pickle.load(f)
f.close()
#print(TS_data)
#print(time_points)
#print(gene_names)
#(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(TS_data, time_points,compute_quality_scores=True,save_models=True)
#print(np.mean(VIM))
exp_data = np.loadtxt('arabidopsis.meristem.expression.csv', delimiter=",", skiprows=1, usecols=list(range(1,21)))
to_remove = []
print(exp_data.shape)
for x in range(exp_data.shape[0]):
    if (np.sum(exp_data[x,(0,2,4,6,8,10,12,14,16,18)])==0) or (np.sum(exp_data[x,(1,3,5,7,9,11,13,15,17,19)])==0):
        to_remove.append(x)
#print(to_remove)
#np.delete(exp_data)
exp_data = np.delete(exp_data,obj=to_remove,axis=0)
print(exp_data.shape)
nulls = 0
for x in exp_data:
    if np.sum(x)==0:
        nulls+=1
print(nulls)
exp1 = exp_data[:,(0,2,4,6,8,10,12,14,16,18)].T
exp2 = exp_data[:,(1,3,5,7,9,11,13,15,17,19)].T
trans_exp_data = [exp1,exp2]
arab_time_points = [np.array(range(1,11)) for x in (1,2)]
#print(arab_time_points)
(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(trans_exp_data, arab_time_points, compute_quality_scores=True)
pd.DataFrame(VIM).to_csv('VIM.csv')
pd.DataFrame(alphas).to_csv('alphas.csv')
#print(alphas)
print(prediction_score)
print(stability_score)
