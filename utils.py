import math 
import numpy as np
def abs_error(actual, correct):
	return abs(actual - correct)

def rel_error(actual, correct):
	return abs_error(actual, correct) / correct

def per_error(actual, correct):
	return rel_error(actual, correct) * 100.0

def euc_distance(x, y):
	return abs(x) * np.sqrt(1.0 + ((y)/(x))**2)
