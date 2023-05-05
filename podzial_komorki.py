import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from decimal import Decimal

n0=1
n1=3.165 - 1j*0.1
n2=3.165 + 1j*0.1
lamb=1.55e-6
k=2*np.pi/lamb


J10 = [[(n1 + n0 )/(2*n1), (n1 - n0)/(2*n1)],
       [(n1 - n0)/(2*n1), (n1 + n0)/(2*n1)]]

J21 =  [[(n2 + n1 )/(2*n2), (n2 - n1)/(2*n2)],
       [(n2 - n1)/(2*n2), (n2 + n1)/(2*n2)]]

J12 = [[(n1 + n2 )/(2*n1), (n1 - n2)/(2*n1)],
       [(n1 - n2)/(2*n1), (n1 + n2)/(2*n1)]]

J02 = [[(n0 + n2 )/(2*n0), (n0 - n2)/(2*n0)],
       [(n0 - n2)/(2*n0), (n0 + n2)/(2*n0)]]

N=26



def C_N1(I_out):
    wynik = [np.sqrt(I_out),0]
    return wynik

lamb_lmab=0.473933
Isg= 1e4
Isa= 1e4
a = lamb_lmab*lamb/2
b = lamb_lmab*lamb/2


def aN_bN_Q(I_out):
    return np.dot(J10,C_N1(I_out))

def cN_dN_Q(aN_bN_0):
    return np.dot(J21,aN_bN_0)


def P1_i(aN_bN,Q):
    ni=0.1
    an_abs_po = np.power(np.absolute(aN_bN[0]), 2)
    bn_abs_po = np.power(np.absolute(aN_bN[1]), 2)
    in1l = ni/(1 + (an_abs_po + bn_abs_po) / Isg)
    P = [
        [ np.exp(-1j * k * ( n1.real - in1l ) * a/Q) ,              0 ],
        [ 0 ,               np.exp(1j * k * (n1.real - in1l) * a / Q)]
    ]
    return P
def P2_i(cN_dN,Q):
    ni=0.1
    cn_abs_po = np.power(np.absolute(cN_dN[0]), 2)
    dn_abs_po = np.power(np.absolute(cN_dN[1]), 2)
    in2l = ni/(1 + (cn_abs_po + dn_abs_po) / Isa)

    P = [
        [ np.exp(-1j * k * ( n1.real + in2l ) * a/Q) ,              0 ],
        [ 0 ,               np.exp(1j * k * (n1.real + in2l) * b / Q)]
    ]

    return P

def I_in_ref(I_out,Q):
    aN_bN = aN_bN_Q(I_out)
    for j in range(1,N):
        for i in range(Q,0,-1):
            aN_bN = np.dot(P1_i(aN_bN,Q),aN_bN)
        cN_dN = cN_dN_Q(aN_bN) # przejście z warstwy

        for i in range(Q,0,-1):
            cN_dN = np.dot(P2_i(cN_dN,Q),cN_dN)
        aN_bN = np.dot(J12,cN_dN)

    a0_b0 = np.dot(J02,cN_dN)

    I_in = np.power(np.absolute(a0_b0[0]), 2)
    I_ref = np.power(np.absolute(a0_b0[1]), 2)
    return I_in, I_ref

def I_in(I_out,Q):
    return I_in_ref(I_out,Q)[0]

def I_ref(I_out,Q):
    return I_in_ref(I_out,Q)[1]

Q = 10 #ilosc podwarstw
I_out_tab = []

for i in range (0,8):
    I_out_temp=np.arange(10**(i),10**(i+1),10**(i-1))
    I_out_tab = np.concatenate((I_out_tab, I_out_temp))


I_in_tab = [I_in(i,Q) for i in I_out_tab]
I_ref_tab = [I_ref(i,Q) for i in I_out_tab]

plt.plot(I_out_tab, I_in_tab)

plt.ylabel(r'$P_{out}[W/m^2]$')
plt.xlabel(r'$P_{in}[W/m^2]$')
plt.yscale('log')
plt.xscale('log')
plt.show()


I_ref_to_I_inc = []
for r,i in zip(I_ref_tab,I_in_tab):
    I_ref_to_I_inc.append(r/i)

plt.plot(I_ref_to_I_inc,I_in_tab)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P_{in}[W/m^2]$')
plt.ylabel("R")
plt.show()

I_out_to_I_inc = []
for o,i in zip(I_out_tab,I_in_tab):
    I_out_to_I_inc.append(o/i)

plt.plot(I_out_to_I_inc,I_in_tab)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P_{in}[W/m^2]$')
plt.ylabel("T")
plt.show()