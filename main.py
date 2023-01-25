import numpy as np
import math
delta_lambda_tab = []

def T(delta_lambda):
    n0=1
    n2=3.165 + 1j*0.1
    n1=3.165 - 1j*0.1

    # Współczynniki r i t przy wejściu z lewej strony
    r02 = (n0-n2)/(n0+n2)
    r20 = (n2-n0)/(n0+n2)
    t20 = 2*n2/(n2+n0)
    t02 = 2*n0/(n0+n2)

    r12 = (n1-n2)/(n1+n2)

    # Współczynniki r i t przy przejściu prz środek komórki PT
    r21 = (n2-n1)/(n2+n1)

    t12 = 2*n1/(n1+n2)

    t21 = 2*n2/(n2+n1)

    r12 = (n1-n2)/(n1+n2)

    #
    r10 = (n1 - n0) / (n1 + n0)

    t01 = 2 * n0 / (n1 + n0)

    t10 = 2 * n1 / (n0 + n1)

    r01 = (n0 - n1) / (n0 + n1)




    S02 = [[r02, t20],
           [t02, r20]]

    S21 = [[r21, t12],
           [t21, r12]]

    S10 = [[r10, t01],
           [t10, r01]]

    t2p=math.exp()


    return np.dot(S10,S21)


print(T(10))