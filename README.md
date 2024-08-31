# two_states_noise_induced_transitions
Répertoir des codes rédigés pour un stage auprès de V.Lucarini à propos des transitions induites par du bruit entre deux états stables compétitifs de modèles bistables inspirés du climat.

**PARTIE 1 : Modèle bistable simple défini à partir d'un potentiel avec deux minimums**

**1.bistable_model** Résolution numérique d'une équation stochastique définie à partir d'un double puits de potentiel et étude du modèle avec calcul du temps moyen passé dans un état stable. 

Utilisation de la méthode de résolution des EDS ITO-Euler avec une bibliothèque adaptée. 
Tracer de trajectoire du système. 
Calcul du nombre de transitions entre états compétitifs et du temps moyen passé dans un état. 
Fit du temps moyen afin de retrouver des valeurs caractéristiques du système

Dans le fichier *"bistable_model"* on procède à la résolution de l'équation différentielle stochastique suivante : $\dot{X} = -4X^3+2\alpha X + \sigma dW$ définie à partir du potentiel $V : X \longmapsto X^4 - \alpha X^2$. Ce potentiel possède deux minimums que l'on calcule grâce à la fonction *argrelextrema* de la bibliothèque SCIPY. Pour résoudre l'équation on fait appel à la fonction *itoint* de la bibliothèque SDEINT qui utilise la méthode de résolution Ito-Euler. Itoint prend comme argument entre autre deux fonctions correspondant à la partie déterministe et à la partie stochastique de l'équation, c'est pourquoi on définit la fonction $X \longmapsto -4X^3+2\alpha X$ plus tôt dans le code. Cette résolution permet déjà de tracer la trajectoire du système $X(t)$. On procède pour cela à un paramétrage important des données de temps d'intégrations et d'intensité du bruit afin d'obtenir une trajectoire lisible et avec de nombreux sauts entre états stables. On peut ensuite tracer la densité de probabilité de présence du système en fonction de la position grâce à la fonction *hist* de la bibliothèque MATPLOTLIB.PYPLOT.

On définit de plus une fonction *count_transition* pour compter à partir d'une trajectoire calculée le nombre de transitions que connaît le système entre les deux états stables. Pour ce faire on utilise un masque sur les données afin de réaliser l'approximation prévue : si le système est proche de l'état stable il est considéré dans l'état stable. Ce masque supprime les valeurs à l'intérieur d'une bande entre les deux états stables et de cette manière supprime les variations parasites qui ne font pas passer d'un état à un autre. Suite à cela grâce à la fonction *groupby* de la bibliothèe ITERTOOLS on regroupe les valeurs entre elles selon la condition : $x>0$. Cette fonction groupe les valeurs à proximité et nous donne donc une séquence de True et de False dont la longueur nous donne le nombre de transitions. 
Finalement on calcule $<\tau>$ et on plot nos valeurs. On réalise un fit de la courbe. 

**2. bistable_model_dv** Résolution numérique d'une équation stochastique définie à partir d'un double puits de potentiel et étude du modèle avec calcul du temps moyen passé dans un état stable. 

Utilisation de la méthode de résolution des EDS ITO-Euler avec une bibliothèque adaptée. 
Tracer de trajectoire du système. 
Calcul du nombre de transitions entre états compétitifs et du temps moyen passé dans un état. 
Fit du temps moyen afin de retrouver des valeurs caractéristiques du système.

Le fichier *"bistable_model_dv"* reprend exactement la même démarche que le précédent à l'exception cette fois que c'est la profondeur du double puits de potentiel qui varie. Rien dans le code des fonctions de comptage et de tracer ne diffère à l'exception que l'on y calcule $\DeltaV$ lorsqu'on fait varier $\alpha$ du potentiel dans la fonction *calcul_taumoy*. 

**PARTIE 2 : Two box model de Stommel bistable**

**1. stommel_model_diagramme** Tracer du diagramme de bifurcation pour le modèle de Stommel étudié. 

Résolution avec la bibliothèque sympy d'équation différentielles déterministes afin de trouver les points fixes du systèmes
Trie des solutions données par la bibliothèque pour ne conserver que les réelles (utilise une tolérance pour ne conserver que les solutions avec une partie complexe très petite)
Tracer du diagramme

Dans le fichier *"stommel_model_diagramme* on procède au tracer du diagramme de bifurcation. Ici la résolution numérique n'est nécessaire que pour éviter un calcul litteral trop lourd. On résoud le système d'équations différencielles déterministes grâce à la fonction *solve* de la bibliothèque SYMPY. Nous ne cherchons que des solutions réelles et cette fonction ne renvoie que des solutions complexes, dont la pluspart n'ont en réalité qu'une partie complexe presque égale à 0. On procède donc au tri de ces solutions grâce à une tolérance que l'on fixe à $10^{-7}$: si la partie complexe de la solution est inférieure à cette tolérance on la conserve. On peut alors tracer le diagramme. 

**2. stommel_integrate** Résolution numérique d'un système d'équations différentielles stochastiques défini à partir du modèle climatique de Stommel et étude du modèle avec calcul du temps moyen passé dans un état stable. 

Utilisation de la méthode de résolution des EDS ITO-Euler avec la bibliothèque sdeint pour résoudre le système de dimension 2
Tracer de trajectoire du système dans l'espace des phases
Tracer de trajectoire du système en 1D
Calcul du nombre de transitions entre états compétitifs et du temps moyen passé dans un état. 
Fit du temps moyen afin de retrouver des valeurs caractéristiques du système.

Le fichier *"stommel_integrate"* comporte la résolution du système d'équations différentielles stochastiques, le tracer des trajectoires et le calcul de $<\tau>$ pour différentes valeurs de $\sigma$ l'intensité du bruit. Pour résoudre l'EDS on fait à nouveau appel à la fonction *itoint* de la bibliothèque SDEINT, mais cette fois on définit les fonctions sur l'espace des phases. La fonction *stommel_uncoupled* défini les équations en une matrice ligne et la fonction *noise* ajoute le bruit (s qu'elle prend en argument) seulement sur l'équation $dS$. On trace ensuite la trajectoire dans l'espace des phases et en 1D. 
On remarque grâce à l'histogramme que les deux états stables n'ont pas le même poids, notre manière de calculer $<\tau>$ change. On compte le temps passé dans l'état de départ. Pour cela dans la fonction *count_transitions* on sélectionne dans les valeurs de $X$ après le masque celles qui sont dans l'états que l'on souhaite avec la condition : $x[mask]>x1+coeffmask$ et l'on compte le nombre de points qui respectent cette condition. Ainsi en le multipliant au pas d'intégration on obtient le temps passé en tout dans cet état. On compte les transtions de la même manière pour autant. On réalise ensuite un fit exponentiel.
