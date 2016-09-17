#!/usr/bin/python
#*- coding: utf-8 -*-

# AUTHOR: Mathieu Tortuyaux
# DATE: 20/06/16
# DERNIERE MAJ: 20/06/16

from __future__ import division
import matplotlib.pyplot as plt
from time import time


# -- DECLARATIONS FONCTIONS 

def feval(funcName, *args):

	#
	# Similaire à `feval` in Matlab
	#
	# Exemple : feval('cos', pi) = -1
	#

	return eval(funcName)(*args)

def volterra(time, y):

	#
	# -- FONCTION IDENTIQUE A CELLE EN MATLAB
	#


	# -- CONSTANTE(S)

	alpha = 0.01


	# -- CALCUL DES DERIVEES
	dy = []

	dy.append(2 * y[0] - alpha * y[0] * y[1])
	dy.append(-y[1] + alpha * y[0] * y[1])  

	return dy

def add(tab, value):

	#
	# -- ADDITIONNE UNE VALEUR SUR TOUS LES ÉLÉMENTS DU VECTEUR
	#

	for i in range(0, len(tab)):


		tab[i] = tab[i] + value[i]

	return tab

def mul(tab, value):

	#
	# -- MULTIPLIE UNE VALEUR SUR TOUS LES ÉLÉMENTS DU VECTEUR
	#


	for i in range(0, len(tab)):

		tab[i] = tab[i] * value

	return tab


def rk4(y, dt, t_final, fonction):


	# -- DECLARATION VARIABLES LOCALES

	time = 0					# START
	n = round(t_final / dt)		# NOMBRE D'ITERATIONS
	t = []						# CREATION DU TABLEAU DE TEMPS
	dy = [0] * 2				# VECTEUR DE SORTIE

	# -- CONDITIONS INITIALES

	t.append(time)
	dy[0] = y[0]
	dy[1] = y[1]

	inter = []					# STOCKER L'EVOLUTION DES RENARDS
	inter1 = []					# STOCKER L'EVOLUTION DES LAPINS

	# -- ITERATIONS DE RUNGE KUTTA

	for i in range(0, int(n)):

		k1 = feval(fonction, time, dy)								#k1 = f(t, y)
		
		k2 = feval(fonction, time+(dt / 2), add(dy, mul(k1, dt/2)))	#k2 = f(t + dt/2, y + k1 * (dt/2))
		
		k3 = feval(fonction, time+(dt / 2), add(dy, mul(k2, dt/2)))	#k3 = f(t + dt/2, y + k2 * (dt/2))

		k4 = feval(fonction, time+dt, add(dy, mul(k3, dt)))			#k4 = f(t + dt, y+k3 )
		
		dy[0] = dy[0] + (dt/6) * (k1[0] + 2*(k2[0]+k3[0]) + k4[0])
		dy[1] = dy[1] + (dt/6) * (k1[1] + 2*(k2[1]+k3[1]) + k4[1])


		time = time + dt
		t.append(time)	



		inter.append(dy[0])
		inter1.append(dy[1])

	t.pop()	# on retire le dernier élément pour que tous les vecteurs aient la même taille


	return t, inter, inter1

if __name__ == "__main__":

	ci = [300, 125]

	start = time()
	vecTemps, vecRenard, vecLapin= rk4(ci, 0.001, 25, 'volterra')
	stop = time()

	print("Calcul réalisé en %s s" %(stop - start))


	plt.plot(vecTemps, vecRenard, vecTemps, vecLapin)
	plt.show()
