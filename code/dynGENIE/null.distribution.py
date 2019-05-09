import numpy as np
import pandas as pd
from dynGENIE3 import *

# Want to make reproducible results
np.random.seed(0)

# First need to import and clean up the data
exp_data = np.loadtxt('arabidopsis.meristem.expression.csv', delimiter=",", skiprows=1, usecols=list(range(1,21)))
print(exp_data.shape)

# Then need to remove the null rows
to_remove = []
for x in range(exp_data.shape[0]):
    if (np.sum(exp_data[x,(0,2,4,6,8,10,12,14,16,18)])==0) or (np.sum(exp_data[x,(1,3,5,7,9,11,13,15,17,19)])==0):
        to_remove.append(x)
exp_data = np.delete(exp_data,obj=to_remove,axis=0)

# Completely Random
cr_exp_data = np.random.permutation(exp_data.flatten()).reshape(exp_data.shape)
# Random Within Rows
rwr_exp_data = np.copy(exp_data)
for x in rwr_exp_data:
    np.random.shuffle(x)
#rwr_exp_data = np.random.permutation(exp_data.T).T

# Again, need to ensure that there is no rows of zeros which will crash the program
cr_to_remove = []
rwr_to_remove = []
for x in range(exp_data.shape[0]):
    if (np.sum(cr_exp_data[x,(0,2,4,6,8,10,12,14,16,18)])==0) or (np.sum(cr_exp_data[x,(1,3,5,7,9,11,13,15,17,19)])==0):
        cr_to_remove.append(x)
        #print("cr:",x)
    if (np.sum(rwr_exp_data[x,(0,2,4,6,8,10,12,14,16,18)])==0) or (np.sum(rwr_exp_data[x,(1,3,5,7,9,11,13,15,17,19)])==0):
        rwr_to_remove.append(x)
        #print("rwr:",x)
cr_exp_data = np.delete(cr_exp_data,obj=cr_to_remove,axis=0)
rwr_exp_data = np.delete(rwr_exp_data,obj=rwr_to_remove,axis=0)

# Want to know if the shapes are the same
print(cr_exp_data.shape,rwr_exp_data.shape)

# Now run the actual analysis
arab_time_points = [np.array(range(1,11)) for x in (1,2)]

exp1 = cr_exp_data[:,(0,2,4,6,8,10,12,14,16,18)].T
exp2 = cr_exp_data[:,(1,3,5,7,9,11,13,15,17,19)].T
trans_cr_exp_data = [exp1,exp2]
(cr_VIM, cr_alphas, cr_prediction_score, cr_stability_score, cr_treeEstimators) = dynGENIE3(trans_cr_exp_data, arab_time_points,compute_quality_scores=True)

pd.DataFrame(cr_VIM).to_csv('cr_VIM.csv')
pd.DataFrame(cr_alphas).to_csv('cr_alphas.csv')

exp1 = rwr_exp_data[:,(0,2,4,6,8,10,12,14,16,18)].T
exp2 = rwr_exp_data[:,(1,3,5,7,9,11,13,15,17,19)].T
trans_rwr_exp_data = [exp1,exp2]
(rwr_VIM, rwr_alphas, rwr_prediction_score, rwr_stability_score, rwr_treeEstimators) = dynGENIE3(trans_rwr_exp_data, arab_time_points,compute_quality_scores=True)

pd.DataFrame(rwr_VIM).to_csv('rwr_VIM.csv')
pd.DataFrame(rwr_alphas).to_csv('rwr_alphas.csv')

# Print out the data, saving that which is needed
print("Completely random:")
print(cr_prediction_score)
print(cr_stability_score)

print("Random within rows:")
print(rwr_prediction_score)
print(rwr_stability_score)
