import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

#X=[1,2,3]
#Z=[1,2,3]
#plt.plot(X, Z)

#print(str(np.arange(0, 2, step=2/kroki)))

def adjust_ticks(max_value_new,max_value_old,steps):

    Ticks=[]
    Tick_Label=[]

    for i in range (0,steps):
        Ticks.append(i*max_value_old/steps)
        Tick_Label.append(i*max_value_new/steps)


    return Tick_List,Tick_Label



def adjust_ticks(max_value_new,steps):
    Tick_Labels=[]
    for i in range (0,steps):
        Tick_Labels.append(i*max_value_new/(steps))

    return Tick_Labels

ax = plt.figure().add_subplot(projection='3d')

# Plot a sin curve using the x and y axes.
x = np.linspace(0, 5, 100)
y = np.sin(x * 2 * np.pi) / 2 + 0.5
ax.plot(x, y, zs=0, zdir='z', label='curve in (x, y)')



labels=adjust_ticks(2,5,10)
ax.set_xticks([0.0,0.2,0.4,0.6,0.8, 1.0,3,4,4.3,5])
ax.set_xticklabels(labels)
plt.show()
