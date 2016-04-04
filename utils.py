import math 
import numpy as np

def abs_error(actual, correct):
	return abs(actual - correct)

def rel_error(actual, correct):
	return abs_error(actual, correct) / correct

def per_error(actual, correct):
	return rel_error(actual, correct) * 100.0

def number_of_floats_in_system(base, sigdigits, m, M):
	return 2(base - 1)(base ^ (sigdigits -1))(M - m + 1) + 1

def euclidean_distance_from_origin(x, y):
	return abs(x) * np.sqrt(1.0 + ((y)/(x))**2)

# Numerically stable quadratic formula
def quadratic_formula(a, b, c):
	x1 = (-b + np.sqrt(b*b - (4 * a * c)))/(2 * a)
	x2 = (-b - np.sqrt(b*b - (4 * a * c)))/(2 * a)

	print(x1)

	if abs(x1) > abs(x2):
		x2 = c / (a * x1)
		return [x1, x2]
	else:
		x1 = c / (a * x2)
		return [x1, x2]

# converts the mantissa to binary equivalent
# m is the number of bits the mantissa occupies
def dec_to_bin_mantissa(mantissa, m):
	r = mantissa
	b = ""
	for k in range(1, m):
		if r >= 2**k:
			b += '1'
			r = r - 2 ** (-1 * k)
		else:
			b += '0'

	return b