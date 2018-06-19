#coding:utf-8


import numpy as np
import scipy.linalg as la
from matplotlib import pyplot as plt

def spline(UG, OG, stellenzahl):
	X = np.linspace(UG, OG, stellenzahl)
	H = [X[i+1]-X[i] for i in range(stellenzahl-1)]
	lamda = [H[i]/(H[i]+H[i+1]) for i in range(0, stellenzahl-2)]
	my = [H[i+1]/(H[i]+H[i+1]) for i in range(0, stellenzahl-2)]
	F1 = [1/(1+x**2) for x in X]
	F2 = [(F1[i]-F1[i+1])/(X[i]-X[i+1]) for i in range(stellenzahl-1)]
	F3 = [(F2[i]-F2[i+1])/(X[i]-X[i+2]) for i in range(stellenzahl-2)]
	M = (1./3.)*np.eye(stellenzahl-2)+(1./6.)*my[1]*np.eye(stellenzahl-2, k=1)+(1./6.)*lamda[1]*np.eye(stellenzahl-2, k=-1)
	vect = la.solve(M, F3)
	vect =np.append(vect, [0])
	vect = np.append([0], vect)
	C = [(F1[i]-F1[i-1])/(H[i])-(1./6.)*H[i]*(vect[i-1]-vect[i]) for i in range(stellenzahl-1)]
	D = [F1[i-1]-(1./6.)*vect[i-1]*H[i]**2 for i in range(stellenzahl-1)]
	C=np.append([0], C)
	D=np.append([0], D)
	print("X ist")
	print(X)
	print("H ist")
	print(H)
	print("lamda ist")
	print(lamda)
	print("my ist")
	print(my)
	print("F1 ist")
	print(F1)
	print("F2 ist")
	print(F2)
	print("F3 ist")
	print(F3)
	print(M)
	print("m0 bis m(n) sind")
	print(vect)
	print("C ist")
	print(C)
	print("D ist")
	print(D)

	return X, H, vect, C, D, stellenzahl



def evaluate(y, X, H, vect, C, D, stellenzahl):
	i=0
	while i < stellenzahl-1:
		if X[i] <= y and X[i+1] >=y:
			break
		else:
			i=i+1
	sy= (1./6.)*vect[i]*H[i]*((y-X[i-1])**3)/H[i]-(1./6.)*vect[i-1]*H[i]*((X[i]-y)**3)/H[i]+C[i]*(y-X[i-1])+D[i]
	return sy


def main():
	print("="*60)
	print("Guten Morgen!")
	print("Implementierung der Erstellung eines kubischen Spline für Funktion f(x)=1/(1+x)^2 im Intervall [-5,5] an 9 äquidistanten Stützstellen.")
	X, H, vect, C, D, stellenzahl = spline(-5,5, 50)
	sol0 = evaluate(0.1, X, H, vect, C, D, stellenzahl)
	print("sol0 ist")
	print(sol0)

	spread=np.linspace(-5,5, 5000)
	evalfn = [1/(1+x**2) for x in spread]
	evalspln = [evaluate(x, X, H, vect, C, D, stellenzahl) for x in spread]
	plt.title("Graphen")
	plt.plot(spread, evalfn, 'ko', label="Funktionsauswertung")
	plt.plot(spread, evalspln, 'bo', label="Spline")
	plt.legend()
	plt.grid()
	plt.savefig("plots/plot.pdf")

if __name__=="__main__":
	main()