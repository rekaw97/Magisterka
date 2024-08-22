
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d
from decimal import Decimal

n0=1
n2=3.165 + 1j*0.1
n1=3.165 - 1j*0.1
lamb=1.55e-6
k=2*np.pi/lamb

# Współczynniki r i t przy wejściu z lewej strony
r02= (n0- n2) / (n2 + n0)
r20 = (n2-n0)/(n0+n2)
t02 = 2*n0/(n0+n2)
t20 = 2*n2/(n2+n0)

r12 = (n1-n2)/(n1+n2)

# Współczynniki r i t przy przejściu prz środek komórki PT
r21 = (n2-n1)/(n2+n1)

t12 = 2*n1/(n1+n2)

t21 = 2*n2/(n2+n1)

r12 = (n1-n2)/(n1+n2)

#Współczynniki r i t przy wejściu z prawej strony
r10 = (n1 - n0) / (n1 + n0)

t01 = 2 * n0 / (n1 + n0)

t10 = 2 * n1 / (n1 + n0)

r01 = (n0 - n1) / (n0 + n1)




S02 = [[r02, t20],
       [t02, r20]]

S21 = [[r21, t12],
       [t21, r12]]

S10 = [[r10, t01],
       [t10, r01]]


def t2p(delta_lambda):
    return np.exp(-1j * k * n2 * (delta_lambda * lamb) / 2)


def t1p(delta_lambda):
    return np.exp(-1j * k * n1 * (delta_lambda * lamb) / 2)

def t2l(delta_lambda):
    return np.exp(1j * k * n2 * (delta_lambda * lamb) / 2)

def t1l(delta_lambda):
    return np.exp(1j * k * n1 * (delta_lambda * lamb) / 2)

J02 = [[(1 / t02), -(r20 / t02)],
       [(r02 / t02), (t20 - (r02 * r20) / t02)]]

J21 = [[(1 / t21), -(r12 / t21)],
       [(r21 / t21), (t12 - (r21 * r12) / t21)]]

J10 = [[(1 / t10), -(r01 / t10)],
       [(r10 / t10), (t01 - (r01 * r10) / t10)]]

J12 = [[(1 / t12), -(r21 / t12)],
       [(r12 / t12), (t21 - (r21 * r12) / t12)]]

def P2(delta_lambda):
    return [[t2p(delta_lambda), 0],
            [ 0, t2l(delta_lambda)]]

def P1(delta_lambda):
    return [[t1p(delta_lambda), 0],
            [ 0, t1l(delta_lambda)]]
def M(delta_lambda,ilosc_komorek):
    lista = [J02]
    for i in range(0,ilosc_komorek):
        lista.extend([P2(delta_lambda),J21,P1(delta_lambda),J12])

    lista.append(P2(delta_lambda))
    lista.append(J21)
    lista.append(P1(delta_lambda))
    lista.append(J10)

    M = [[1, 0],
         [ 0, 1]]

    for i in lista:
        M = np.dot(M,i)

    return M
def S(delta_lambda,ilosc_komorek):
    M_ = M(delta_lambda,ilosc_komorek)
    S = [ [( M_[1][0] / M_[0][0]), (M_[1][1]-M_[0][1]*M_[1][0]/M_[0][0])],
                [ (1/M_[0][0]),-(M_[0][1]/M_[0][0])]]

    return S

def RL(delta_lambda,ilosc_komorek):
    return np.power(np.absolute(S(delta_lambda,ilosc_komorek)[0][0]),2)

def TL(delta_lambda,ilosc_komorek):
    return np.power(np.absolute(S(delta_lambda,ilosc_komorek)[1][0]), 2)
def TR(delta_lambda,ilosc_komorek):
    return np.power(np.absolute(S(delta_lambda,ilosc_komorek)[0][1]), 2)
def RR(delta_lambda,ilosc_komorek):
    return np.power(np.absolute(S(delta_lambda,ilosc_komorek)[1][1]), 2)

def rL(delta_lambda,ilosc_komorek):
    return S(delta_lambda,ilosc_komorek)[0][0]

def tL(delta_lambda,ilosc_komorek):
    return S(delta_lambda,ilosc_komorek)[1][0]

def tR(delta_lambda,ilosc_komorek):
    return S(delta_lambda,ilosc_komorek)[0][1]

def rR(delta_lambda,ilosc_komorek):
    return S(delta_lambda,ilosc_komorek)[1][1]


def wart_wlas_1(delta_lambda,ilosc_komorek):
    return 0.5*(rL(delta_lambda,ilosc_komorek)+rR(delta_lambda,ilosc_komorek)+np.sqrt(
                                                            np.power(rL(delta_lambda,ilosc_komorek),2)
                                                            + np.power(rR(delta_lambda,ilosc_komorek),2)
                                                            + 4 * np.power(tR(delta_lambda,ilosc_komorek),2)
                                                            - 2 * rL(delta_lambda,ilosc_komorek) * rR(delta_lambda,ilosc_komorek)))

def wart_wlas_2(delta_lambda,ilosc_komorek):
    return 0.5*(rL(delta_lambda,ilosc_komorek)+rR(delta_lambda,ilosc_komorek)-np.sqrt(
                                                            np.power(rL(delta_lambda,ilosc_komorek),2)
                                                            + np.power(rR(delta_lambda,ilosc_komorek),2)
                                                            + 4 * np.power(tR(delta_lambda,ilosc_komorek),2)
                                                            - 2 * rL(delta_lambda,ilosc_komorek) * rR(delta_lambda,ilosc_komorek)))


delta_lambda_tab = np.arange(0,20,0.01)
RL_tab = []
TL_tab = []
TR_tab = []
RR_tab = []
WW_tab = []
WW1_tab = []
WW2_tab = []
WW_jedn = []



for i in delta_lambda_tab:
    RL_tab.append(RL(i,1))
    TL_tab.append(TL(i,1))
    TR_tab.append(TR(i,1))
    RR_tab.append(RR(i,1))
    WW1_tab.append(np.absolute(wart_wlas_1(i,1)))
    WW2_tab.append(np.absolute(wart_wlas_2(i,1)))


for i in WW_tab:
    if i<1:
        WW1_tab.append(i)
    if i==1:
        WW1_tab.append(i)
        WW2_tab.append(i)
    if i>1:
        WW2_tab.append(i)


for i in range(0,len(delta_lambda_tab)):
    if WW1_tab[i]>WW2_tab[i]:
        WW1_tab[i], WW2_tab[i] = WW2_tab[i], WW1_tab[i]


plt.plot(delta_lambda_tab[1:],TR_tab[1:],label=r'$T_{R}$')
plt.plot(delta_lambda_tab[1:],TL_tab[1:],label=r'$T_{L}$')
plt.yscale('log')
plt.ylabel("T")
plt.xlabel("Λ/λ")
plt.title("Transmitancja")
plt.legend()
plt.show()

plt.plot(delta_lambda_tab[1:],RR_tab[1:],label=r'$R_{R}$')
plt.plot(delta_lambda_tab[1:],RL_tab[1:],label=r'$R_{L}$')
plt.yscale('log')
plt.ylabel("R")
plt.xlabel("Λ/λ")
plt.title("Reflektancja")
plt.legend()
plt.show()


plt.plot(delta_lambda_tab,WW1_tab,label=r'|$λ_{1}$|')
plt.plot(delta_lambda_tab,WW2_tab,label=r'|$λ_{2}$|')
plt.yscale('log')
plt.ylabel("λ")
plt.xlabel("Λ/λ")
plt.title("Wartości własne macierzy S(Λλ)")
plt.legend()
plt.show()


#%%
delta_lambda_tab = np.arange(0,2,0.0001)


licz_kom=150
maksimum =0
maksima_tab=[]
Z_tab=[]
for i,dl in enumerate(delta_lambda_tab):
    Z=[]
    max_temp=0
    for j in range(1,licz_kom+1):
        val = TL(dl, j)
        Z.append(val)


    max_temp=max(Z)
    if max_temp > maksimum:
        maksimum = max_temp
        max_del_lam = dl
    maksima_tab.append(max_temp)
    Z_tab.append(Z)

#%%
def adjust_ticks(max_value_new,max_value_old,steps):

    Ticks_List=[]
    Tick_Label=[]

    for i in range (0,steps+1):
        Ticks_List.append(i * max_value_old / steps)
        Tick_Label.append(i * max_value_new / steps)


    return Ticks_List,Tick_Label[::-1]


kolory = ['green','yellow','orange','red']
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

alphas = maksima_tab/maksimum
alphas_log=np.log10(alphas)
alphas_norm = alphas_log-min(alphas_log)
alphas_line = alphas_norm/max(alphas_norm)
#print(alphas_line)

for i in range(0,len(Z_tab)):


    z=Z_tab[i]

    X=np.arange(1,licz_kom+1,1)
    kol_num = int(np.floor(alphas_line[i])*(len(kolory)-1))
    #print(kol_num)
    if np.max(z)>4:
    #kol_num=int(np.ceil(maksima_tab[i]/maksimum*len(kolory)))-1


        kolor = kolory[kol_num]
        ax.plot(X, z, len(Z_tab)-i, zdir='y', color=kolor,alpha = 1 )


    else:
        kolor = kolory[kol_num]
        ax.plot(X, z, len(Z_tab)-i, zdir='y', color = 'blue',alpha=0.01)

delta_lambda_end=2
ticks,labels=adjust_ticks(delta_lambda_end,len(delta_lambda_tab),5)
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
plt.ylabel("Λ/λ")
plt.xlabel("N")
plt.title(r"$T_L$"+' dla  '+r"$n_1={}$".format(n1))
print(max_del_lam)
print(maksimum)

plt.show()




