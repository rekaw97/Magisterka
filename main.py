import numpy as np
import math
delta_lambda_tab = []

def T(delta_lambda):
    n0=1
    n2=3.165 + 1j*0.1
    n1=3.165 - 1j*0.1
    lamb=1.55*10^(-6)
    k=2*math.pi/lamb

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

    t2p = math.exp(-1j * k * n2 * (delta_lambda * lamb) / 2)
    t1p = math.exp(-1j * k * n1 * (delta_lambda * lamb) / 2)
    t2l = math.exp(1j * k * n2 * (delta_lambda * lamb) / 2)
    t1l = math.exp(1j * k * n1 * (delta_lambda * lamb) / 2)

    J02 = [[(1 / t02), -(r20 / t02)],
           [(r02 / t02), (t20 - (r02 * r20) / t02)]]

    J21 = [[(1 / t21), -(r12 / t21)],
           [(r21 / t21), (t12 - (r21 * r12) / t21)]]

    J10 = [[(1 / t10), -(r01 / t10)],
           [(r10 / t10), (t01 - (r01 * r10) / t10)]]

    P2 = [[t2p, 0],
           [ 0, t2l]]

    P1 = [[t1p, 0],
           [ 0, t1l]]

    lista = [J02,J21,J10,P2,P1]

    M = [[1, 0],
           [ 0, 1]]

    for i in lista:
        M = np.dot(M,i)

    S= [ ( [M[1,0] / M[0,0]]), (M[1,1]-M[0,1]*M[1,0]/M[0,0]),
           [ (1/M[0,0]),-(M[0,1]/M[0,0])]]
    return np.dot(S10,S21)


print(T(10))