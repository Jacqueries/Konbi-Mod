""" Fichier contenant la classe Weights qui définit les poids associés aux modes à optimiser 
"""

import numpy as np
import sys, os
import random
import copy
class Weights:
	"""Classe des poids. Permet de ponderer les modes de basses fréquence
	lors de l'optimisation.
	"""
	def __init__(self,nmodes):
		self.weights = np.full(nmodes,0.0) #initWeights
		self.linf = np.full(len(self.weights),-1.0)
		self.lsup = np.full(len(self.weights),1.0)
		self.wp = copy.deepcopy(self.weights)
		self.combinaisonMax = []
		self.importance = np.full(nmodes,False)
		self.timeCollect = []
		self.memory = self.initMemory()

	def initWeightsO(self,nmodes):
		"""Initialise un vecteur de taille nmodes 
		revoie le vecteur des poids (aleatoirement tiré entre -1 et 1)
		"""
		w = []
		for i in range(nmodes):
			w.append(random.uniform(-1.0,1.0))#*(1-(i/nmodes)))
		return w

	def watchmode(self):
		indice = 0
		for i in range(len(self.weights)):
			if indice != i:
				self.weights[i] = 0.
			else:
				self.weights[i] = 1

	def ajustIndiv(self,indice):
		"""Au début du processus iteratif de Monte Carlo on initialise chque poid 
		par 1 et -1 et les autres par 0.0 pour evaluer la contribution du mode
		"""
		if indice % 2 == 1:
			w = -1
		else:
			w = 1
		self.weights[int(indice/2)] = w

	def reajustWeights(self,indice):
		"""Tire au hasard les poids à l'interieur de leurs bornes respectives
		"""
		if indice < len(self.weights)*2:
			self.weights = np.full(len(self.weights),0.0)
			self.ajustIndiv(indice)
		else:
			self.wp = copy.deepcopy(self.weights)
			for i in range(len(self.weights)):
				if self.importance[i]:
					self.weights[i] = random.uniform(self.linf[i],self.lsup[i])
			self.updateMemory()

	def precState(self):
		"""Return to last save
		"""
		self.weights = copy.deepcopy(self.wp)

	def reajustLimits(self,indice):
		"""Reajuste les limites de tirage aleatoire des poids
		Si un poids a augmenté et que le volume du canal aussi alors sa limite 
		inferieure devient sa valeur précedente et inversement
		"""
		if indice < len(self.weights)*2:
			self.importance[int(indice/2)] = True
		else:
			for i in range(len(self.weights)):
				if self.weights[i] > self.wp[i]:
					self.linf[i] = self.weights[i]
				elif self.weights[i] < self.wp[i]:
					self.lsup[i] = self.weights[i]
	def rdz(self):
		"""Réinitialise un interval choisi aleatoirement
		"""
		i = random.randint(0,len(self.weights)-1)
		self.linf[i] = -1.0
		self.lsup[i] = 1.0


	def saveCombinaisonMax(self,time):
		"""Savegarde les poids qui maximisent le volume d'une conformation
		"""
		self.combinaisonMax.append(copy.deepcopy(self.weights))
		self.timeCollect.append(time)

	def initMemory(self):
		"""Recuperation de données
		"""
		memi = []
		for _ in range(len(self.weights)):
			memi.append([])
		return memi

	def updateMemory(self):
		"""Sauve ma valeur des poids a chaque itération
		"""
		for i in range(len(self.memory)):
			self.memory[i].append(self.weights[i])


def watch(Vol,Svol,w):
	"""Verifie que les valeurs récentes du volume/surface ne stagnent pas
	Si elles stagnent, perturbe n poids tirés aleatoirement ( réinitialise les bornes)
	les valeurs 10 et 0.1 sont arbitraires
	"""
	Svol.append(Vol)
	if len(Svol)>10: # si stagnation pendant 10 itérations
		if np.std(Svol) < 0.1*np.mean(Svol): # si les valeurs varient peu
			n = len(w.weights) # on réinitialise n poids aléatoires
			for i in range(n):
				w.rdz()
		Svol = []
	return [Svol,w]


