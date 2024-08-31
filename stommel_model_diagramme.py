import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from sympy import symbols, Eq, solve, S 

#données
sigma = np.linspace(0.5, 1.2, 30) #interval de sigma considéré

tolerance = 10**(-7) #tolérence pour la sélection des solutions considérées réelles données par le module sympy
psi = []
sig = []

mu = 7.5
theta = 1
ep_a = 0.34

#recherche des points fixes
for s in sigma : 
    #résolution système d'équation pour T - S > 0
    T, S = symbols('T S')
    eq1 = Eq(-(T-theta)/ep_a - T - mu*(T - S)*T, 0)
    eq2 = Eq(s - S - mu*(T - S)*S, 0)
    solution_pos = solve((eq1, eq2), (T, S))
    #résolution système d'équation pour T - S < 0 
    T, S = symbols('T S')
    eq1 = Eq(-(T-theta)/ep_a - T - mu*(S - T)*T, 0)
    eq2 = Eq(s - S - mu*(S - T)*S, 0)
    solution_neg = solve((eq1, eq2), (T, S))
    #suppression des parties complexes très petite dans les solutions données par le module sympy
    solution_pos = [(complex(x).real, complex(y).real) for (x, y) in solution_pos if all(abs(complex(val).imag) < tolerance for val in (x, y))]
    psi_pos = [x - y for (x, y) in solution_pos if x - y >0]
    solution_neg = [(complex(x).real, complex(y).real) for (x, y) in solution_neg if all(abs(complex(val).imag) < tolerance for val in (x, y))]
    psi_neg = [x - y for (x, y) in solution_neg if x - y <0]
    res = psi_pos + psi_neg
    #récupération des résultats en s'assurant que les tableaux aient la même dimension
    for k in range(len(res)):
        sig.append(s)
        psi.append(res[k])
    


#tracer du diagramme
plt.scatter(sig, psi, label = r'$\psi(\sigma) = S - T$ (avec $\mu =$' + str(mu) + r', $\epsilon_a =$' + str(ep_a) + ')')
plt.legend()
plt.ylabel(r'$\psi$')
plt.xlabel(r'$\sigma$')
plt.axis()
plt.title('Diagramme de Bifurcation avec une région de bistabilité')
plt.show()