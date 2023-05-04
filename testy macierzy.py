import numpy as np
tab1 = [[5,2],[3,4]]

tab2 = [3,2]

tab3 = np.transpose([3,2])

print(np.dot(tab1,tab2))
print(np.dot(tab1,tab3))