import numpy as np
import linearSystems as ls
def Vandermonde(X):
	n = X.shape[0]
	A = np.vstack([np.power(X, i).reshape((1, n)) for i in range(n)])
	return A.T


def Monomial(points, x=None):
	A = Vandermonde(points[:, 0])
	y = points[:, 1]
	C = ls.GE(A, y)
	
	if x is None:
		return C

	return np.sum(C[i]*(x**i) for i in range(len(C)))

def L(X, x, k):
	return np.product(np.array([(x - X[i])/(X[k] - X[i])  for i in range(len(X)) if i != k]))

def Lagrangian(points, x):
	X  = points[:, 0]
	Fx = points[:, 1]

	return np.sum(Fx[i] * L(X, x, i) for i in range(len(X)))

def w(X, x):
	return np.product(x - X)

def a(fx, px, wx):
	return (fx - px) / wx

def Newton(points):
	n = len(points)
	X = points[:, 0]
	Fx = points[:, 1]
	coeffs = np.zeros(n, dtype=float)

	coeffs[0] = Fx[0]

	for i in range(1, n):
		fx = Fx[i]
		px = NewtonEval(coeffs[:i], X[:i], X[i])
		wx = w(X[:i], X[i])
		newCoeff = a(fx, px, wx)
		coeffs[i] = newCoeff

	return coeffs

def NewtonEval(coeffs, X, x):
	n = len(coeffs)
	W = np.ones(n, dtype=float)

	for i in range(1, n):
		W[i] = w(X[:i], x)

	return coeffs.dot(W)