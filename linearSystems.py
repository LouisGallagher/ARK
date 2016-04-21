import utils
import numpy as np

def BackSubstitution(A, b):
	n = A.shape[0]
	x = np.zeros(n, dtype=float)

	x[n-1] = b[n-1] / A[n-1, n-1]

	for i in range(n-1, -1, -1): 
		x[i] = (b[i] - sum(A[i, i+1:n] * x[i+1:])) / A[i,i]

	return x

def ForwardSubstitution(A, b):
	n = A.shape[0]
	x = np.zeros(n, dtype=float)

	x[0] = b[0] / A[0, 0]

	for i in range(1, n): 
		x[i] = (b[i] - sum(A[i, 0:i] * x[0:i])) / A[i,i]

	return x

def GE(A, b):
	n = A.shape[0]
	
	for i in range(0, n):
		j = np.argmax(abs(A[i:, i])) + i
		A[[i,j], :] = A[[j, i], :]
		b[i], b[j]= b[j] ,b[i]

		for k in range(i+1, n):
			coeff = A[k, i]/A[i, i]
			
			A[k, i:] = A[k, i:] - (A[i, i:] * coeff)
			b[k] = b[k] - (b[i] * coeff)
			
	return BackSubstitution(A, b)

def LU(A, b = None):
	n = A.shape[0]
	P = np.identity(n)
	L = np.identity(n, dtype=float)
	U = np.zeros((n, n), dtype=float)
	for i in range(0, n):
		j = np.argmax(abs(A[i:, i])) + i
		A[[i,j], :] = A[[j, i], :]
		P[[i, j], :] = P[[j,i], :]
		
		for k in range(i+1, n):
			coeff = A[k, i]/A[i, i]
			
			A[k, i] = A[k, i]/A[i, i]
			A[k, i+1:] = A[k, i+1:] - (A[i, i+1:] * A[k, i])
			
	L += np.tril(A, -1)
	U += np.triu(A)

	if b is None:
		return (U, L, P)

	y = ForwardSubstitution(L, P.dot(b.transpose()).transpose())

	x = BackSubstitution(U, y)
	
	return (U, L, P, x)


def Jacobi(A, b, X0, delta):
	n = A.shape[0]
	XCurr = np.copy(X0)
	XNew = np.zeros(n, dtype=float)
	
	while True:
	 	for i in range(n):
	 		XNew[i] = (b[i] - np.sum(XCurr * A[i, :]) + (XCurr[i] * A[i, i]))/A[i,i]
		
	 	err = np.linalg.norm(XCurr - XNew) / (np.linalg.norm(XNew) + np.finfo(float).eps)
	 	print err
	 	if err < delta:
	 		break

	 	XNew, XCurr = XCurr, XNew

	return XCurr

def GaussSeidel(A, b, X0, delta):
	n = A.shape[0]
	XCurr = np.copy(X0)
	XNew = np.zeros(n, dtype=float)

	while True:
		for i in range(n):
			XNew[i] = (b[i] - np.sum(XNew[0:i] * A[i, 0:i]) - np.sum(XCurr[i+1:] * A[i, i+1:])) / A[i, i]

		err = np.linalg.norm(XCurr - XNew) / (np.linalg.norm(XNew) + np.finfo(float).eps)
		
	 	if err < delta:
	 		break

	 	XNew, XCurr = XCurr, XNew

	return XCurr

def GaussSeidelWithRelaxation(A, b, X0, delta, k=10, p=1):
 	n = A.shape[0]
 	XOld = np.zeros(n, dtype=float)
 	XCurr = np.copy(X0)
 	XNew = np.zeros(n, dtype=float)

 	for i in range(k):
 		for j in range(n):
 			XNew[j] = (b[j] - np.sum(XNew[0:j] * A[j, 0:j]) - np.sum(XCurr[j+1:] * A[j, j+1:])) / A[j, j]

 	 	XOld = XCurr
 	  	XCurr, XNew = XNew, XCurr

 	deltaXKIterations = np.linalg.norm(XOld - XCurr)

 	for i in range(p):
 		for j in range(n):
 			XNew[j] = (b[j] - np.sum(XNew[0:j] * A[j, 0:j]) - np.sum(XCurr[j+1:] * A[j, j+1:])) / A[j, j]
	 	
 	 	XOld = XCurr
 	  	XCurr, XNew = XNew, XCurr

 	deltaXKPlusPIterations = np.linalg.norm(XOld - XCurr)
 	omegaOpt = 2.0 / (1.0 + np.sqrt(1.0 - np.power(deltaXKPlusPIterations / deltaXKIterations, 1.0/p))) 

 	print omegaOpt
 	it = 11
 	while True:
 		for i in range(n):
 			XNew[i] = (b[i] - np.sum(XNew[0:i] * A[i, 0:i]) - np.sum(XCurr[i+1:] * A[i, i+1:])) / A[i, i]
 			XNew[i] = omegaOpt * XNew[i] + (1 - omegaOpt)*XCurr[i]

 		err = np.linalg.norm(XCurr - XNew) / (np.linalg.norm(XNew) + np.finfo(float).eps)
		#print err
		it += 1
 	 	if err < delta:
	 		break

	 	XCurr, XNew = XNew, XCurr
	print it
	return XCurr

def HilbertMatrix(n):
	M = np.zeros((n,n), dtype=float)

	for i in range(n):
		for j in range(n):
			M[i, j] = 1.0 / (float(i) + float(j) + 1.0) 

	return M