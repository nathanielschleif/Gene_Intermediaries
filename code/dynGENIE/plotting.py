import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

test = pd.read_csv("test_VIM.csv",index_col=0,header=0,dtype=float)
print(test.iloc[0])
flt = test.values.flatten()
#flt = flt[flt>0.02]
print(np.amax(test.values.flatten(),axis=None))
#plt = ser.plot.hist(grid=True, bins=10)
plt.hist(flt,bins=40)
plt.title('Histogram of Trimmed dynGENIE3 Scores')
plt.xlabel('Score')
plt.ylabel('Count')
#plt.grid(axis='y', alpha=0.75)
plt.savefig('testdata_histogram.png')
