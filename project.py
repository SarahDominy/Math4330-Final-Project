import BasicFunction.py as BF

def vandermonde(x):
	"""
	Computes the vandermonde matrix out of a vector.

	This algorithm takes in a vector v and computes a vandermonde matrix with a degree of 4 by creating a temporary
	vector, which holds the computations of all the degrees in a value of x, and then adds the temporary vector into
	the matrix A.

	Args:
	   x: A list of numbers, representing a vector.
	Returns:
	   A: A list of lists of numbers, representing the vandermonde matrix created from vector x. The lists are stored
		as the rows of the matrix.
	"""
	A = []
	for iterator in range(len(x)):
		temp = []
		temp.append(x[iterator])
		for i in range(5):
			temp.append(x[iterator]**i)
		A.append(temp)
	return A

def modGS(A):
	"""
	Computes the modified gram-schmidt of a vandermonde matrix.

	This algorithm obtains the QR factorization of a vandermonde matrix A by using a modified version of the
	Gram-Schmidt algorithm. Note: See Take Home exam 1 for the two-norm and orthoDecomp functions.

	Args:
		A: A list of lists of numbers, representing the vandermonde matrix A.
	Returns:
		A matrix that is the result of applying the Gram-Schmidt algorithm to matrix A.
	"""
	for k=1 in range(n):
		v[k] = A[k]
	for k=1 in range(n):
		r[k][k] = BF.twoNorm(v[k])
		q[k] = v[k] * (1/r[k][k])
		for j = k+1 in range(n):
			r[k][j] = BF.dot(BF.conjugateTranspose(q[k]), v[j])
			v[j] = v[j] - r[k][j]*q[k]

def dataFitting(x, y):
	"""
	Generates a set of numbers that can solve the equation Ax = y, where A is the set of numbers to solve for.

	This algorithm creates a set of numbers that, if multiplied by the corresponding elements in vector x, will
	result in a number equal to the corresponding element in vector y. This is done by first generating a vandermonde
	matrix from vector x, then applying the modified Gram-Schmidt algorithm to the matrix to obtain the QR factorization.
	Then, the algorithm computes the inverse of Q, Q*, and uses backsubstitution to solve for Ra = Q*y, where a is unknown,
	but will be used to obtain the degree for interpolating polynomial.

	Args:
                x: A list of numbers, representing a matrix that contains the numbers for the input values of a function.
		y: A list of numbers, representing a matrix that contains the numbers for the output values of a function.
	Returns:
		A list of numbers, representing the numbers whose corresponding values will, when multiplied with the
			corresponding values in matrix x, will result in the corresponding values in matrix y.
	"""
	A = vandermonde(x)
	modGS(A)
	Q = BF.inverse(Q)
	alpha = BF.backSub(Q, y)
	return alpha
