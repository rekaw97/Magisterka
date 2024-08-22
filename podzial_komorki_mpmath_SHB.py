import numpy as np
import mpmath as mp
from matplotlib import pyplot as plt
from matplotlib.ticker import NullFormatter
import mpl_toolkits.axisartist as AA

mp.mp.dps=20
Im = 0.1
n0=1
n1 = mp.mpc( 3.165, -Im)
n2 = mp.mpc(3.165 , Im)

lamb=0.000155
k=mp.fdiv(2*np.pi,lamb)



#lamb_lmab=0.473933
lamb_lmab=1.42047
#lamb_lmab=0.7892





a = lamb_lmab*lamb/2
b = lamb_lmab*lamb/2
alfa = 0

def sat_coeff(wektor_amp,Is,SH):


    SHB =SH *( wektor_amp[0] * mp.conj(wektor_amp[1]) +  mp.conj(wektor_amp[0]) * wektor_amp[1])

    result = mp.fdiv(1,
                     1 +
                        ((mp.power(mp.fabs(wektor_amp[0]), 2) + mp.power(mp.fabs(wektor_amp[1]), 2)
                                        + SHB ) / Is))

    return result


def j10(n1u,C_N1,ni,SH):


    J10 = mp.matrix( [[(n1u + n0) / (2 * n1u), (n1u - n0) / (2 * n1u)],
            [(n1u - n0) / (2 * n1u), (n1u + n0) / (2 * n1u)]])
    aN_bNi =  J10 * C_N1

    in1l = mp.re(ni * sat_coeff(aN_bNi, Isg, SH))


    n1u = mp.mpc(3.165, - in1l)

    return aN_bNi , n1u

def j20(n2u,A_N1,ni,SH):


    J20 = mp.matrix( [[(n2u + n0) / (2 * n2u), (n2u - n0) / (2 * n2u)],
            [(n2u - n0) / (2 * n2u), (n2u + n0) / (2 * n2u)]])
    cN_dNi =  J20 * A_N1

    in1l = mp.re(ni * sat_coeff(cN_dNi, Isa, SH))


    n2u = mp.mpc(3.165, in1l)

    return cN_dNi , n2u

#macierz przejścia z ośrodka 1 do ośrodka 2, dynamiczna na podstawie c i d
def j21(n1u,n2u,aN_bN,ni,SH):

    J21 = mp.matrix([[(n2u + n1u) / (2 * n2u), (n2u - n1u) / (2 * n2u)],
            [(n2u - n1u) / (2 * n2u), (n2u + n1u) / (2 * n2u)]])

    cN_dNi = J21 * aN_bN


    in2l = mp.re(ni * sat_coeff(cN_dNi, Isa, SH))

    n2u = mp.mpc(3.165 , in2l)


    return cN_dNi, n2u

def j12(n1u,n2u,cN_dN,ni,SH):

    J12 =mp.matrix( [[(n1u + n2u) / (2 * n1u), (n1u - n2u) / (2 * n1u)],
        [(n1u - n2u) / (2 * n1u), (n1u + n2u) / (2 * n1u)]])

    aN_bNi = J12 * cN_dN

    in1l = mp.re(ni * sat_coeff(aN_bNi, Isg, SH))

    n1u = mp.mpc(3.165, - in1l)


    return aN_bNi, n1u

def j02(n2u):
    J02 = mp.matrix([[(n0 + n2u )/(2*n0), (n0 - n2u)/(2*n0)],
            [(n0 - n2u)/(2*n0), (n0 + n2u)/(2*n0)]])

    return J02


def j01(n1u):
    J01 = mp.matrix([[(n0 + n1u )/(2*n0), (n0 - n1u)/(2*n0)],
            [(n0 - n1u)/(2*n0), (n0 + n1u)/(2*n0)]])

    return J01

def C_N1(I_out):
    wynik = mp.matrix([mp.sqrt(I_out),0])
    return wynik


def P1_i(aN_bN,Q,n1u,ni,SH):

    P = mp.matrix( [
            [ mp.exp(-1j * k * n1u * a/Q) ,              0 ],
            [ 0 ,               mp.exp(1j * k * n1u * a / Q)]
        ])

    aN_bNi = P * aN_bN

    n1l = mp.re(ni * sat_coeff(aN_bNi, Isg, SH))

    n1u = mp.mpc(3.165, - n1l)

    return aN_bNi,n1u
def P2_i(cN_dN,Q,n2u,ni,SH):

    P =mp.matrix( [
            [ mp.exp(-1j * k * n2u * b/Q) ,              0 ],
            [ 0 ,               mp.exp(1j * k * n2u * b/Q)]
         ])

    cN_dNi = P * cN_dN

    n2l = mp.re(ni * sat_coeff(cN_dNi, Isa, SH))

    n2u = mp.mpc(3.165 , n2l)

    return cN_dNi,n2u




N=21
def I_in_ref_L(I_out,Q,n1_,n2_,SH): #główna funkcja zwracająca natężenie fali padającej oraz odbitej

    n1u = n1_
    n2u = n2_
    ni = Im

    for g in range(0, 4):  # uzgodnienie współcznynnika załamania w warstwie 1
        aN_bN,n1u = j10(n1u,C_N1(I_out),ni,SH)


    for j in range(0,N): #pętla po komórkach

        for i in range(Q,0,-1): #pętla po podwarstwach w warstwie 1
            aN_bN,n1u=P1_i(aN_bN, Q, n1u,ni,SH)

        n2u = n2_


        for g in range(0,4): #uzgodnienie współczynnika załamania w warstwie 2
            cN_dN,n2u = j21(n1u,n2u,aN_bN,ni,SH)


        for i in range(Q,0,-1): #pętla po podwarstwach w warstwie 2
            cN_dN,n2u = P2_i(cN_dN, Q, n2u,ni,SH)

        n1u = n1_
        for g in range(0,4): #uzgodnienie współczynnika załamania w warstwie 1
            aN_bN,n1u = j12(n1u,n2u,cN_dN,ni,SH)


    J02 = j02(n2u)
    a0_b0 = J02 * cN_dN #przejście z warstwy 2 do powietrza

    I_in = mp.power(mp.fabs(a0_b0[0]), 2)
    I_ref = mp.power(mp.fabs(a0_b0[1]), 2)
    return I_in, I_ref



def I_in_ref_R(I_out,Q,n1_,n2_,SH): #główna funkcja zwracająca natężenie fali padającej oraz
    n1u = n1_
    n2u = n2_
    ni = Im

    for g in range(0, 4):  # uzgodnienie współcznynnika załamania w warstwie 1
        cN_dN,n2u = j20(n2u,C_N1(I_out),ni,SH)


    for j in range(0,N): #pętla po komórkach

        for i in range(Q,0,-1): #pętla po podwarstwach w warstwie 1
            cN_dN,n2u=P2_i(cN_dN, Q, n2u,ni,SH)

        n1u = n1_


        for g in range(0,4): #uzgodnienie współczynnika załamania w warstwie 2
            aN_bN,n1u = j12(n1u,n2u,cN_dN,ni,SH)


        for i in range(Q,0,-1): #pętla po podwarstwach w warstwie 2
            aN_bN,n1u = P1_i(aN_bN, Q, n1u,ni,SH)

        n2u = n2_
        for g in range(0,4): #uzgodnienie współczynnika załamania w warstwie 1
            cN_dN,n2u = j21(n1u,n2u,aN_bN,ni,SH)


    J01 = j01(n1u)
    c0_d0 = J01 * aN_bN #przejście z warstwy 2 do powietrza

    I_in = mp.power(mp.fabs(c0_d0[0]), 2)
    I_ref = mp.power(mp.fabs(c0_d0[1]), 2)
    return I_in, I_ref



parametry = []
parametry.append([1e2,1e2,1e-3,1e6])
parametry.append([1e2,1e4,1e-4,8*1e7])
parametry.append([1e2,1e6,1e-4,1e10])
parametry.append([1e4,1e2,1e-4,1e8])
parametry.append([1e4,1e4,1e-1,8*1e7])
parametry.append([1e4,1e6,1e-2,1e10])
parametry.append([1e6,1e2,1e-4,1e10])
parametry.append([1e6,1e4,1e-2,8*1e9])
parametry.append([1e6,1e6,1e0,1e10])
for par in parametry:
    Isg= par[0]
    Isa= par[1]

    start=par[2]
    stop= par[3]


    Q = 10 #ilosc podwarstw

    I_out_tab = []

    x=mp.power(mp.fdiv(stop, start), mp.fdiv(1, 2000))
    for i in range(0,2000):#wartości natężenia na wyjściu
        I_out_tab.append(start * mp.power(x,i))


    I_in_tab_L = []
    I_ref_tab_L =[]
    I_in_tab_L_SHB = []
    I_ref_tab_L_SHB =[]

    I_in_tab_R = []
    I_ref_tab_R =[]
    I_in_tab_R_SHB = []
    I_ref_tab_R_SHB =[]

    for i in I_out_tab:
        I_in_L,I_ref_L=I_in_ref_L(i, Q, n1, n2, 0)
        I_in_R,I_ref_R=I_in_ref_R(i, Q, n2, n1, 0)
        I_in_L_SHB,I_ref_L_SHB=I_in_ref_L(i, Q, n1, n2, 1)
        I_in_R_SHB,I_ref_R_SHB= I_in_ref_R(i, Q, n2, n1, 1)


        I_in_tab_L.append(I_in_L)
        I_ref_tab_L.append(I_ref_L)
        I_in_tab_R.append(I_in_R)
        I_ref_tab_R.append(I_ref_R)

        I_in_tab_L_SHB.append(I_in_L_SHB)
        I_ref_tab_L_SHB.append(I_ref_L_SHB)
        I_in_tab_R_SHB.append(I_in_R_SHB)
        I_ref_tab_R_SHB.append(I_ref_R_SHB)


    ##%%
    plt.figure(figsize=(11, 5))
    plt.plot(I_in_tab_L, I_out_tab,label='$I_{out}^{(a)}$ η=0',linestyle='-',linewidth=2,color='navy')
    plt.plot(I_in_tab_R, I_out_tab,label='$I_{out}^{(g)}$ η=0',linestyle='-',linewidth=2,color='orange')
    plt.plot(I_in_tab_L_SHB, I_out_tab,label='$I_{out}^{(a)}$ η=1',linestyle=':',linewidth=2,color='limegreen')
    plt.plot(I_in_tab_R_SHB, I_out_tab,label='$I_{out}^{(g)}$ η=1',linestyle=':',linewidth=2,color='red')


    plt.ylabel(r'$I_{out}^{(\{g,a\})}$' +' '+ '$[W/cm^2]$',fontsize="15")
    plt.xlabel(r'$I_{in}$' +' '+ '$[W/cm^2]$',fontsize="15")
    plt.yscale('log')
    plt.xscale('log')
    plt.gca().tick_params(axis='both', which='major', labelsize=15)

    plt.legend(fontsize="15")

    plt.figtext(0.8,0.3,'$I_{sg}$'+' ' +'$=10^{}W/cm^2$'.format( int(np.log10(Isg)))+"\n"+
                '$I_{sa}$'+' ' +'$=10^{}W/cm^2$'.format( int(np.log10(Isa))),fontsize="15",bbox=dict(facecolor='white', alpha=0.8,edgecolor='lightgrey'))

    plt.grid(True, which="both", ls="-")

    plt.gca().xaxis.set_minor_locator(plt.LogLocator(base=10, subs='all', numticks=100))
    plt.gca().axes.get_xaxis().set_minor_formatter(NullFormatter())
    plt.margins(x=0)

    plt.tight_layout()

    plt.savefig('I_out Isg={:.0e} Isa={:.0e}.png'.format(Isg,Isa))
    plt.show()


    T_tab_L= []
    T_tab_L_SHB= []
    for o,i in zip(I_out_tab, I_in_tab_L): #wyliczenie T
        T_tab_L.append(mp.fdiv(o, i))

    for o,i in zip(I_out_tab, I_in_tab_L_SHB): #wyliczenie T
        T_tab_L_SHB.append(mp.fdiv(o, i))


    T_tab_R = []
    T_tab_R_SHB = []
    for o,i in zip(I_out_tab, I_in_tab_R_SHB): #wyliczenie T
        T_tab_R.append(mp.fdiv(o, i))
    for o,i in zip(I_out_tab, I_in_tab_R_SHB): #wyliczenie T
        T_tab_R_SHB.append(mp.fdiv(o, i))


    plt.figure(figsize=(11, 5))
    plt.plot(I_in_tab_L, T_tab_L, label='$T_a$ η=0',linestyle='-',linewidth=2,color='navy')
    plt.plot(I_in_tab_R, T_tab_R, label='$T_g$ η=0',linestyle='-',linewidth=2,color='orange')
    plt.plot(I_in_tab_L_SHB, T_tab_L_SHB, label='$T_a$ η=1',linestyle=':',linewidth=2,color='limegreen')
    plt.plot(I_in_tab_R_SHB, T_tab_R_SHB, label='$T_g$ η=1',linestyle=':',linewidth=2,color='red')

    plt.ylabel(r'$T_a,T_g$' +' '+ '$[j.u.]$',fontsize="15")
    plt.xlabel(r'$I_{in}$' +' '+ '$[W/cm^2]$',fontsize="15")
    plt.yscale('log')
    plt.xscale('log')
    plt.gca().tick_params(axis='both', which='major', labelsize=15)
    plt.legend(fontsize="15")

    plt.figtext(0.15,0.3,'$I_{sg}$'+' ' +'$=10^{}W/cm^2$'.format( int(np.log10(Isg)))+"\n"+
        '$I_{sa}$'+' '+'$=10^{}W/cm^2$'.format( int(np.log10(Isa))),fontsize="15",bbox=dict(facecolor='white', alpha=0.8,edgecolor='lightgrey'))


    plt.grid(True,axis='x', which="both", ls="-")
    plt.grid(True,axis='y', which="major", ls="-")
    plt.gca().xaxis.set_minor_locator(plt.LogLocator(base=10, subs='all', numticks=100))
    plt.gca().axes.get_xaxis().set_minor_formatter(NullFormatter())
    plt.margins(x=0)
    #plt.title('$I_{sg}$'+'=$10^{}$'.format( int(np.log10(Isg)))+'$[W/cm^2]$'+' $I_{sa}$'+'=$10^{}$'.format(int(np.log10(Isa)))+'$[W/cm^2]$')
    plt.tight_layout()
    plt.savefig('T Isg={:.0e} Isa={:.0e}.png'.format(Isg,Isa))
    plt.show()

    ########################################

    R_tab_L = []
    R_tab_L_SHB = []
    for r,i in zip(I_ref_tab_L, I_in_tab_L): #wyliczenie R
        R_tab_L.append(mp.fdiv(r, i))

    for r, i in zip(I_ref_tab_L_SHB, I_in_tab_L_SHB):  # wyliczenie R
         R_tab_L_SHB.append(mp.fdiv(r, i))




    R_tab_R = []
    R_tab_R_SHB = []
    for r,i in zip(I_ref_tab_R, I_in_tab_R): #wyliczenie R
        R_tab_R.append(mp.fdiv(r, i))

    for r,i in zip(I_ref_tab_R_SHB, I_in_tab_R_SHB): #wyliczenie R
        R_tab_R_SHB.append(mp.fdiv(r, i))

    plt.figure(figsize=(11,5))
    plt.plot(I_in_tab_L, R_tab_L, label='$R_a$ η=0',linestyle='-',linewidth=2,color='navy')
    plt.plot(I_in_tab_R, R_tab_R, label='$R_g$ η=0',linestyle='-',linewidth=2,color='orange')
    plt.plot(I_in_tab_L_SHB, R_tab_L_SHB, label='$R_a$ η=1',linestyle=':',linewidth=2,color='limegreen')
    plt.plot(I_in_tab_R_SHB, R_tab_R_SHB, label='$R_g$ η=1',linestyle=':',linewidth=2,color='red')
    plt.yscale('log')
    plt.xscale('log')

    plt.ylabel(r'$R_a,R_g$' +' '+ '$[j.u.]$',fontsize="15")
    plt.xlabel(r'$I_{in}$' +' '+ '$[W/cm^2]$',fontsize="15")
    plt.gca().tick_params(axis='both', which='major', labelsize=15)
    plt.legend(fontsize="15")
    plt.figtext(0.15,0.3,'$I_{sg}$'+' ' +'$=10^{}W/cm^2$'.format( int(np.log10(Isg)))+"\n"+
    '$I_{sa}$'+' ' +'$=10^{}W/cm^2$'.format( int(np.log10(Isa))),fontsize="15",bbox=dict(facecolor='white', alpha=0.8,edgecolor='lightgrey'))

    plt.grid(True,axis='x', which="both", ls="-")
    plt.grid(True,axis='y', which="major", ls="-")
    plt.gca().xaxis.set_minor_locator(plt.LogLocator(base=10, subs='all', numticks=100))
    plt.gca().axes.get_xaxis().set_minor_formatter(NullFormatter())
    plt.margins(x=0)
    #plt.title('$I_{sg}$'+'=$10^{}$'.format( int(np.log10(Isg)))+'$[W/cm^2]$'+' $I_{sa}$'+'=$10^{}$'.format(int(np.log10(Isa)))+'$[W/cm^2]$')
    plt.tight_layout()
    plt.savefig('R Isg={:.0e} Isa={:.0e}.png'.format(Isg,Isa))
    plt.show()
