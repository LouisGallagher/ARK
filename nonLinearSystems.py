import numpy as np
import linearSystems as ls

def Jacobian(fprime, x):
	return np.array([[fprime[j, i](x) for i in range(fprime.shape[1])] for j in range(fprime.shape[0])], dtype=float)


def NewtonRaphson(f, fprime, X0, delta):
	n = len(f)
	X1 = np.zeros(n, dtype=float)
	fx = [f[i](X0) for i in range(n)]
	
	while not(np.allclose(fx, np.zeros(n), delta)):
		deltaX = ls.GE(Jacobian(fprime, X0), -1.0 * X0)
		X1 = X0 + deltaX
		X0 = X1
		fx = [f[i](X0) for i in range(n)]
		print fx

	return X0

# def systemToFunction(A, b):
# 	return lambda x: A.dot(x)-b

# def Jacobian(f, x, h):
# 	return np.array([[(f[i,j](x+h)-f[i,j](x-h))/(2.0 * h) for i in range(fprime.shape[1])] for j in range(fprime.shape[0])], dtype=float)

# def NewtonRaphson(f, X0, delta):
# 	n = len(b)
# 	X1 = np.zeros(n, dtype=float)
# 	fx = np.array(f[i](X0) for i in range(len(x)))

# 	while np.allclose(fx, np.zeros(n), delta):
# 		deltaX = ls.GE(Jacobian(f, X0), -1.0 * X0, 0.01)
# 		X1 = X0 + deltaX
# 		X0 = X1
# 		fx = f(X0)

# 	return X0





