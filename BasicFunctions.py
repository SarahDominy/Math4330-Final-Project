def abso(sca):
    """Calculate the absolute value

    Checks if the input is complex or negative. If it complex it takes the square root of the sum of squares of the real and imaginary part.
    If it is negative it multiplies the number by -1.

    Args:
        sca: a number that can be real or complex

    Returns:
        a number which is the absolute value of th
    """

    absVal = 0
    if (sca.imag != 0):
        absVal = (sca.real ** 2 + sca.imag ** 2) **(1/2)
    elif (sca < 0):
        absVal = - sca
    else
        absVal = sca
    return absVal

def twoNorm(vec):
    """Calculates the 2-norm of a vector

    Sums the square of the elements of a given vector and returns the root of the sum.

    Args:
        vec: a list of real and complex numbers representing a vector

    Returns:
        A scalar which is the 2-norm of the given vector
    """

    norm = 0
    for i in range(len(v)):
        norm = norm + (absVal(v[i]) ** 2)
    norm = norm ** (1/2)
    return norm

def scalVecMulti(sca, vec):
    """Multiples together the scalar and the vector

    Goes element by element and takes the product of the scalar and the number in that element

    Args:
        sca: an int that is a scalar
        vec: a list of real and complex numbers representing a vector

    Returns:
        A vector where each element is the product of the scalar and the correspoding element of the original vector
    """
    proVec=vec
    for i in range(len(vec)):
        proVec[i] = sca * vec[i]
    return proVec

def normalize(vec, n):
    """Normalizes a vector with respect to the n norm

    Checks to see if a norm is 0, and if not, divides the matrix by the nth norm

    Args:
        vec: a list of real and complex numbers representing a vector
        n: an int represents what norm is being requested

    Returns:
        a normalized version of the input vector
    """

    norm = twoNorm(vec)

    if (norm == 0):
        print("Invalid Input")
    elif (norm == 1):
        return vec
    else:
        return scalVecMulti((1/norm), vec)

def dot(vec1, vec2):
    """Takes the dot product of the two vectors

    Goes element by element and take the sum of the products of the elements with the corresponding index

    Args:
        vec1: a list of real and complex numbers representing a vector
        vec2: a list of real and complex numbers representing a vector

    Returns:
        A scalar that is the sum of the products of both vectors
    """

    pro = 0
    for i in range(len(vec1)):
        pro = pro + vec1[i] * vec2[i]
    return pro

def vecAdd(vec1,vec2):
    """Adds together two vectors

    Goes element by element of both matrices and adds them together

    Args:
        vec1: a list of real and complex numbers representing a vector
        vec2: a list of real and complex numbers representing a vector or same length as vec1
    Returns:
        A vector that each element is the sum of the 2 input vectors

    """

    sumv=vec1
    if (len(vec1) != len(vec2)):
        print("Invalid inputs, vectors are not of same length")
    else:
        for i in range(len(vec1)):
            sumv[i] = vec1[i] + vec2[i]
        return sumv

def conjugate(compNum):
    """Compute the conjuate of the complex number.

        This function takes a complex number as its input. Initializes a complex number 0+0j as its output and then replaces the real part of its output with the real part of its input and replaces the imaginary part of its output with the imaginary part of its input.

        Args:
            compNum: a complex number
        Returns: 
            the complex conjugate of our input"""
    if compNum.imag==0:
        return compNum
    else:
        conjCompNum=compNum.real
        conjCompNum=conjCompNum-compNum.imag*1j
        return conjCompNum

def conjugateTranspose(A):
	"""
	Computes the conjugate transpose of a matrix.

	This function takes in a matrix and computes the conjugate transpose by first creating a transpose of matrix A, and then
	replacing each element of the transpose matrix with its conjugate value.

	Args:
		A: A list of lists, which contain numbers, representing a matrix.
	Returns:
		A matrix whose rows are the columns of A, and whose columns are the rows of A. Each element in this matrix should
		also be the conjugate of it's corresponding value in the matrix A.
	"""
	result = transpose(A)
	for iterator in range(len(A))
		for element to (len(A[0]))
                        result[iterator][element] = conjugate(result[iterator][element]
	return result


def transpose(A):
	"""
	Computes the transpose of a matrix.

	This function takes in a matrix and computes its transpose by storing each row of the input matrix into a
	temporary value, then placing the temporary value into the correct position on the result matrix.

	Args:
		A: A list of lists, which contain numbers, representing a matrix.
	Returns:
		A matrix whose rows are the columns of matrix A, and whose columns are the rows of matrix A.
	"""
	result = []
	for iterator in range(len(A[0])
		temp = []
		for element in range(len(A))
			temp.append(A[iterator][element])
		result.append(temp)
	return result

def backSub(A,b):
	"""
	Solves the Ax=b equation by taking a matrix A and vector b and finding the x values.
	
	Solves the last row of A and b, then uses that and backsubs it into the above row, repeating until completely solving.
	
	Args:
		A: a list of lists representing a uppertriangular matrix 
		b: a list representing a vector
	Returns:
		result: a vector containing the x values
	"""
	result = []
	for iterator in range(len(A[0]))
		a = len(A[0])
		k = a – iterator + 1
		for j in range(k , (len(A)):
			e = scalVecMulti(A[a – iterator][k], result[k])
                        y = y + (scalVecMulti((A[a – iterator][a – iterator]**-1), e))
		result [a – iterator] = (b[a – iterator] – y)
	return result


def inverse(Q):
	"""
	Computes the inverse of a matrix Q.

	This algorithm takes in a unitary matrix Q and computes its inverse, Q*, by first replacing the rows of Q with
	the columns of Q, then computing the conjugate of each entry in the matrix.

	Args:
		Q: A list of lists of numbers, representing the unitary matrix Q.
	Returns:
		Q*, which is the inverse of matrix Q.
	"""
	result = Q
	for iterator in range(len(Q)):
		for element in range(len(Q[0])):
			result[iterator][element] = Q[element][iterator]
	for iterator in range(len(Q)):
		for element in range(len(Q[0])):
			result[iterator][element] = result[iterator][element]
	return result
