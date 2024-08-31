import numpy as np
import matplotlib.pyplot as plt
import random as rd
import sdeint as sd
from scipy.signal import argrelextrema
import itertools as it
from scipy.optimize import curve_fit


#minimums potentiel
"""amin = 0
amax = 5
Nb_val_a = 300

a = np.linspace(amin, amax, Nb_val_a)
alpha = 0
beta = 0 #le potentiel est symétrique
"""
def potentiel(x, a) : 
    return x**4-a*x**2

pot = np.frompyfunc(potentiel, 1, 1) #transforme le potentiel en ufunc numpy


def min_pot(alpha):
    """entrée : un potentiel avec deux minimums
    sortie : les coordonnées des minimums du potentiel
    """
    x = np.linspace(-2,2,300)
    pot = np.frompyfunc(lambda x : x**4-alpha*x**2 , 1, 1)
    i1 = argrelextrema(pot(x), np.less)[0][0] #indice min 1
    i2 = argrelextrema(pot(x), np.less)[0][1] #indice min 2
    return x[i1], x[i2]

def deltaV(a) : 
    """entrée : une valeur pour alpha coefficient dans le potentiel
    sortie : deltaV, profondeur du puits pour ce potentiel
    """
    x1, x2 = min_pot(a)
    return potentiel(0,a) - potentiel(x1,a)


#fonctions pour ITO-euler et tau moyen

def f(x, t, alpha) : #opposé de la dérivée du potentiel
    return -4*x**3+2*alpha*x



def euler(sigma, x_0, T, nb_val, alpha):
    #euler
    temps = np.linspace(0, T, nb_val) #intervalle des temps parcourus pendant Euler
    return sd.itoint(lambda x, t : -4*x**3+2*alpha*x , lambda x, t : sigma, x_0, temps, np.random)


def count_transitions(x, coeff_mask, x1, x2):
    """entrée : tableau x de la trajectoire, coefficient de masquage, état stable 1, état stable2
    sortie : le nombre de transitions entre états sur la trajectoire
    """
    mask = np.logical_or(x < x1 + coeff_mask, x > x2 - coeff_mask) #masque pour filtré les valeurs de x
    return len([key for key, gp in it.groupby(x[mask]>0)])

def calcul_taumoy(amn, amidmn, amidmx, amx, nb_val_a1, nb_val_a2, nb_val_a3, t, nb_val_t, coeff_mask, sigma): 
    ''' entrée : bornes de l'intervalle d'alpha, nb de valeurs de alpha à considérer, 
    temps d'intégration, nb de séquences pour l'intervalle de temps, coefficient de masquage,
    intensité du bruit
    sortie : tableau de valeurs de deltaV, tableau de valeurs de taumoy correpondantes
    appel à : count_transitions, min_pot, potentiel, euler 
    '''
    alpha1 = np.linspace(amn, amidmn, nb_val_a1, endpoint = False)[1:]
    alpha2 = np.linspace(amidmn, amidmx, nb_val_a2, endpoint = False)
    alpha3 = np.linspace(amidmx, amx, nb_val_a3) 
    alpha = np.concatenate((alpha1, alpha2, alpha3)) #intervalle des sigmas parcourus pour tracer tau
    taumoy = []
    dv = []
    for a in alpha : 
        x_1, x_2 = min_pot(a)
        taumoy.append(T/count_transitions(euler(sigma, x_1, t, nb_val_t, a), coeff_mask, x_1, x_2))
        dv.append(- potentiel(x_1, a)) #c'est bien deltaV car il n'y a pas de terme cst donc V(0) = 0
    return np.array(dv), np.array(taumoy)



#données
sigma = 0.6
#pour euler
T = 15000
Nb_val_t = 1000000


#intervalle de alpha considéré
amin = 0
amidmin = 2*(0.3)**(1/2)
amidmax = 2*(0.7)**(1/2)
amax = 2*(0.9)**(1/2)
Nb_val_a1 = 5
Nb_val_a2 = 10
Nb_val_a3 = 3

ep = 0.1 #coeff de sélection pour le masquage dans count_transitions



dv, taumoy = calcul_taumoy(amin, amidmin, amidmax, amax, Nb_val_a1, Nb_val_a2, Nb_val_a3, T, Nb_val_t, ep, sigma)
plt.scatter(dv, taumoy)
plt.plot(dv, taumoy, label = 'sigma = ' + str(sigma)+ ', période pour euler =' +str(T)+', pas de Euler ='+str(T/Nb_val_t)+', nombre de points=' +str(Nb_val_a1+Nb_val_a2+Nb_val_a3))
plt.xlabel("Delta V")
plt.ylabel("tau moyen")
plt.title('Tau moyen en fonction de deltaV')
plt.legend()
plt.show()

#fit exponentiel

def expo_fit(x, a, b):
    return a*np.exp(b*x)

params, _ = curve_fit(expo_fit, dv, taumoy, p0=[0, 0.1])

dV_fit = np.linspace(0, 0.9, 300)
tau_fit = expo_fit(dV_fit, params[0], params[1])


plt.scatter(dv, taumoy, label = r'$<\tau>(\alpha)$ (avec $\sigma = 0.6$)')
plt.plot(dV_fit, tau_fit, color = 'red', label = r'fit exponentiel  $a\cdot e^{b \cdot x}$ : a=' + str(params[0]) + ', b=' + str(params[1]))
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$<\tau>$')
plt.title("Temps moyen passé dans un des états stables en fonction de la profondeur des puits de potentiel")
plt.legend()
plt.axis()
plt.show()