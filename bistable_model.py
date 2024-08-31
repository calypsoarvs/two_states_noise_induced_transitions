import numpy as np
import matplotlib.pyplot as plt
import random as rd
import sdeint as sd
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
import itertools as it


#minimums potentiel
def min_pot(f) : 
    """entrée : un potentiel avec deux minimums
    sortie : les coordonnées des minimums du potentiel
    """
    x = np.linspace(-2,2,300)
    i1 = argrelextrema(f(x), np.less)[0][0] #indice min 1
    i2 = argrelextrema(f(x), np.less)[0][1] #indice min 2
    return x[i1], x[i2]


#fonctions pour ITO-euler et tau moyen

def f(x, t) : #opposé de la dérivée du potentiel
    return -4*x**3+2*alpha*x-beta



def count_transitions(x, coeff_mask, x1, x2):
    """entrée : tableau x de la trajectoire, coefficient de masquage, état stable 1, état stable2
    sortie : le nombre de transitions entre états sur la trajectoire
    """
    mask = np.logical_or(x < x1 + coeff_mask, x > x2 - coeff_mask) #masque pour filtré les valeurs de x
    return len([key for key, gp in it.groupby(x[mask]>0)])

def calcul_taumoy(sigmn, sigmidmn, sigmidmx, sigmx, nb_val_s1, nb_val_s2, nb_val_s3, temps, f, x_0, coeff_mask, x1, x2): 
    """entrée : bornes de l'intervalle de sigma, nombre de valeurs pour les intervalles de sigma, 
                tableau temps pour l'intégration, fonction dans l'équation diff, 
                état de départ pour l'intégration, coefficient de masquage, états stables 1 et 2
        sortie : tableau des valeurs de sigma, tableau des valeurs de <tau> correspondantes """
    sig1 = np.linspace(sigmn, sigmidmn, nb_val_s1, endpoint = False) 
    sig2 = np.linspace(sigmidmn, sigmidmx, nb_val_s2, endpoint = False)
    sig3 = np.linspace(sigmidmx, sigmx, nb_val_s3) 
    sig = np.concatenate((sig1, sig2, sig3)) #intervalle des sigmas parcourus pour tracer tau
    taumoy = [T/count_transitions(sd.itoint(f, lambda x, t : s, x_0, temps, np.random), coeff_mask, x1, x2) for s in sig]
    return sig, np.array(taumoy)



#potentiel
alpha = 1
beta = 0

def potentiel(x) : 
    return x**4-alpha*x**2+beta*x
pot = np.frompyfunc(potentiel, 1, 1) #transforme le potentiel en ufunc numpy


#données


min1, min2 = min_pot(pot) #min du potentiel

#pour ITO-euler
T = 10000
Nb_val_t = 1000000
temps = np.linspace(0, T, Nb_val_t)

#intervalle de sigma considéré
sigmin = 0.15
sigmidmin = 0.2
sigmidmax = 0.6
sigmax = 0.8
Nb_val_s1 = 3
Nb_val_s2 = 10
Nb_val_s3 = 3
ep = 0.1 #coeff de sélection pour le masquage


x_0 = min1 #position de départ pour la résolution ITO-Euler

#tracer d'une trajectoire pour un sigma précis 
sigma = 0.3
x_rd = sd.itoint(f, lambda x, t : sigma, x_0, temps, np.random)
mask = np.logical_or(x_rd < min1+ep, x_rd > min2-ep)
x_filtré = x_rd[mask]
x_filtré = x_filtré>0
nb_et = len([key for key, gp in it.groupby(x_filtré)])
print(nb_et)
tau = T/nb_et
plt.plot(temps, x_rd, label = 'sigma = ' + str(sigma)+ ', pas=' + str(T/Nb_val_t)+', nb de transisition='+str(nb_et))
plt.legend()
plt.title('X en fonction du temps')
plt.show()

#histogramme densité de probabilité pour un tracer de trajectoire

plt.hist(x_rd, bins = 1000, label = 'sigma = ' + str(sigma) + ", pas d'intégration=" + str(T/Nb_val_t))
plt.legend()
plt.xlabel('position')
plt.ylabel('P')
plt.title('Densité de probabilité en fonction de la position')
plt.show()

#tracer de la loi
sig, taumoy = calcul_taumoy(sigmin, sigmidmin, sigmidmax, sigmax, Nb_val_s1, Nb_val_s2, Nb_val_s3, temps, f, x_0, ep, min1, min2)

plt.scatter(sig, taumoy)
plt.plot(sig, taumoy, label = 'période pour euler =' +str(T)+',pas de Euler ='+str(T/Nb_val_t)+',nombre de points=' +str(Nb_val_s1 + Nb_val_s2 + Nb_val_s3))
plt.xlabel("sigma")
plt.ylabel("tau moyen")
plt.title('Tau moyen en fonction de sigma')
plt.legend()
plt.show()

#fit pour la loi 
def expo_fit(x, a, b):
    return a*np.exp(b/x**2)

params, _ = curve_fit(expo_fit, sigma, taumoy, p0=[0.1, -0.01])

sig_fit = np.linspace(0.03, 0.047, 100)
tau_fit = expo_fit(sig_fit, params[0], params[1])

plt.scatter(sigma, taumoy, label = r'$<\tau>(\sigma)$ (avec $\sigma_0 = 0.9$)')
plt.plot(sig_fit, tau_fit, color = 'red', label = r'fit exponentiel décroissante $a\cdot e^{\frac{b}{x^2}}$ : a=' + str(params[0]) + ', b=' + str(params[1]))
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$<\tau>$')
plt.title("Temps moyen passé dans un des états stables en fonction de l'intensité du bruit")
plt.legend()
plt.show()
