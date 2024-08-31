import numpy as np
import sdeint as sd
import matplotlib.pyplot as plt
import itertools as it
from scipy.optimize import curve_fit

#données 

    #données des équations couplées
mu = 7.5
theta = 1
ep_a = 0.34
sigma_0 = 0.9
sigma = 0.05 #intensité du bruit pour le tracer d'une seule trajectoire

    #intervalle de temps pour l'intégration
tmax = 10000
nb_val_T = 1000000 
tspan = np.linspace(0, tmax, nb_val_T) 

    #états stables dans l'espace des phases
u_0 = np.array([0.57, 0.4]) # on state + état initial pour la résolution
u_0_off = np.array([0.72, 0.74]) # off state

    #états stables pour psi
min2 = 0.57 - 0.4 
min1 = 0.72 - 0.74

#fonctions pour ITO-EUler

def stommel_uncoupled(u, t):
    #définition des équations comme une matrice ligne
    T, S = u[0], u[1]
    dT = -(T - theta) / ep_a - T - mu*np.abs(T - S) * T
    dS = sigma_0 - S - mu*np.abs(T - S) * S
    return np.array([dT, dS])

def noise(u, t) :
    #fonction qui ajoute du bruit seulement à la deuxième équation
    return np.array([[0, 0], [0, sigma]])


#résolution pour une seule valeur de sigma
res = sd.itoint(stommel_uncoupled, noise, u_0, tspan, np.random) #résolution des équations avec la méthode ITO-Euler
T = [T for [T, S] in res]
S = [S for [T, S] in res]

#tracer de la trajecoire dans l'espace des phases pour un sigma
plt.plot(S, T, label = r' T(S) (pour $\sigma_0$ = ' + str(sigma_0) + r" et $\sigma$ = " + str(sigma) + ')')
plt.legend()
plt.title("Trajectoire dans l'espace des phases")
plt.xlabel('S')
plt.ylabel('T')
plt.show()

#tracer de la trajectoire pour psi pour un sigma
psi_values = [T - S for [T, S] in sd.itoint(stommel_uncoupled, noise, u_0, tspan, np.random)]

plt.plot(tspan, psi_values, label = r'$\mu$' + ' = ' + str(mu) + r', $\theta$' + ' = ' + str(theta) + r', $\sigma$ = ' + str(sigma) + r', $\sigma_0$ = ' + str(sigma_0))
plt.legend()
plt.xlabel('temps')
plt.ylabel(r'$\psi$')
plt.title(r'$\psi$ en fonction du temps')
plt.axhline(min1, c='red', linewidth = 0.5)
plt.axhline(min2, c='red', linewidth = 0.5)
plt.show()

#histogramme : densité de probabilité de présence du système dans un état psi
plt.hist(psi_values, bins = 1000, label = r'$P(\psi)$')
plt.legend()
plt.xlabel(r'$\psi')
plt.ylabel('P')
plt.title(r'Densité de probabilité en fonction de $\psi$ ')
plt.show()

#fonctions de calcul de taumoy

def count_transitions(x, coeff_mask, x1, x2):
    """entrée : tableau x de la trajectoire, coefficient de masquage, état stable 1, état stable2
    sortie : le nombre de transitions entre états sur la trajectoire
    """
    mask = np.logical_or(x < x1 + coeff_mask, x > x2 - coeff_mask) #masque pour filtré les valeurs de x
    x_select = x[mask]>x1+coeff_mask
    temps_en_haut = (T/Nb_val_t) * list(x_select).count(True) # multiplication du pas d'intégration par le nombre de points 'en haut'
    nb_en_haut = len([key for key, gp in it.groupby(x_select) if key])
    return temps_en_haut/nb_en_haut


def calcul_taumoy(sigmn, sigmidmn, sigmidmx, sigmx, nb_val_s1, nb_val_s2, nb_val_s3, temps, f, x_0, coeff_mask, x1, x2): 
    """entrée : bornes de l'intervalle de sigma, nombre de valeurs pour les intervalles de sigma, 
                tableau temps pour l'intégration, fonction dans l'équation diff, 
                état de départ pour l'intégration, coefficient de masquage, états stables 1 et 2
        sortie : tableau des valeurs de sigma, tableau des valeurs de <tau> correspondantes """
    sig1 = np.linspace(sigmn, sigmidmn, nb_val_s1, endpoint = False) 
    sig2 = np.linspace(sigmidmn, sigmidmx, nb_val_s2, endpoint = False)
    sig3 = np.linspace(sigmidmx, sigmx, nb_val_s3) 
    sig = np.concatenate((sig1, sig2, sig3)) #intervalle des sigmas parcourus pour tracer tau
    taumoy = [count_transitions(np.array([T - S for [T, S] in sd.itoint(f, lambda x, t : np.array([[0, 0], [0, s]]), x_0, temps, np.random)]), coeff_mask, x1, x2) for s in sig]
    return sig, np.array(taumoy)


#tracer de la loi exponentielle 

#pour euler-ITO
T =20000
Nb_val_t = 1000000
temps = np.linspace(0, T, Nb_val_t)

#intervalle de sigma considéré
sigmin = 0.03
sigmidmin = 0.03
sigmidmax = 0.04
sigmax = 0.045
Nb_val_s1 = 0
Nb_val_s2 = 7
Nb_val_s3 = 3
ep = 0.03#coeff de sélection pour le masquage


u_0 = np.array([0.57 , 0.4]) #état de départ


sig, taumoy = calcul_taumoy(sigmin, sigmidmin, sigmidmax, sigmax, Nb_val_s1, Nb_val_s2, Nb_val_s3, temps, stommel_uncoupled, u_0, ep, min1, min2)
plt.scatter(sig, taumoy)
plt.plot(sig, taumoy, label = 'période pour euler =' +str(T)+',pas de Euler ='+str(T/Nb_val_t)+',nombre de points=' +str(Nb_val_s1 + Nb_val_s2 + Nb_val_s3))
plt.xlabel("sigma")
plt.ylabel("tau moyen")
plt.title('Tau moyen en fonction de sigma')
plt.legend()
plt.show()

#fit 
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
