import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd


overall=pd.read_csv('Overall_Ref_Set_1.csv',header=None)
overall.columns = ['Annual', 'MinRev', 'Hedge', 'Fund']

max_complexity=overall[(overall.Hedge)==1]
min_complexity=overall[(overall.Hedge)==0]


#sns.set(font_scale=2)
#sns.set_style("white")
#sns.set_style("white")

sns.set_context("paper", font_scale=2)

sns.distplot(max_complexity.Annual , color="red", kde=False, bins=10)
sns.distplot(min_complexity.Annual , color="skyblue", kde=False, bins=10)
#plt.xlabel('Maximum Reserve Fund Balance', fontsize=14)
#plt.ylabel('Frequency', fontsize=14)

#plt.legend()


