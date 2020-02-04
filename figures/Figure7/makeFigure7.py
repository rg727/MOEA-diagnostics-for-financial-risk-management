import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
 
objects = ('Borg', 'MOEAD', 'NSGAII', 'NSGAIII', 'RVEA')
objects = ('Borg', 'NSGA-II', 'NSGA-III', 'RVEA', 'MOEA/D')
y_pos = np.arange(len(objects))
l = [0.812751923781605, 0.08794430194210333, 0.05643092707951631, 0.008794430194210334,0.06339318431659949]
l = [x * 100 for x in l]
plt.bar(y_pos, l, align='center',color=['Blue', 'Orange', 'Green', 'Purple', 'Red'])
plt.xticks(y_pos, objects)
plt.ylabel('Contribution to Reference Set (%)',fontsize=14)
plt.xlabel('Algorithm',fontsize=16)
plt.tick_params(axis = 'both', which = 'major', labelsize = 14)


plt.title('Percent Contribution of each MOEA to the Reference Set',fontsize=15)
 
plt.show()