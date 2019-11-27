'''
Calculate the hip joint centre, as defined by Gamage, Lasenby J. (2002)

'''

__author__ = "Andrea Cereatti, Alex Woodall"

import numpy as np
from distance import distance

def calc_HJC(TrP):
	'''
	This function calculates the hip joint centre (HJC)

	Input: TrP clean matrix containing markers'trajectories in the proximal system of reference.
           dim(TrP)=Nc*3p where Nc is number of good samples and p is the number of distal markers

	Output: Cm vector with the coordinates of hip joint center (Cx,Cy,Cz).

	Comments: metodo1b extracts HJC position as the centre of the optimal spherical suface that minimizes the root 
			  mean square error between the radius(unknown) and the distance of the centroid of marker's coordinates 
			  from sphere center(unknown). Using definition of vector differentiation is it possible to put the 
			  problem in the form: A*Cm=B that is a linear equation system

	References: Gamage, Lasenby J. (2002). 
            	New least squares solutions for estimating the average centre of rotation and the axis of rotation.
            	Journal of Biomechanics 35, 87-93 2002   
                Halvorsen correzione bias
	'''

	r, c = np.size(TrP)
	D = np.zeros(3)
	V2 = np.array([])
	b1 = np.array([0, 0, 0])

	for j in range (1, 3, c):
		d1 = np.zeros(3)

		for i in range(1, r):
			d1 = np.array([d1 + np.transpose(TrP[i-1,(j-1):(j+1)]) * TrP[i-1,(j-1):(j+1)]])
			a = np.power(TrP[(i-1),(j-1)],2) + np.power(TrP[i-1,j],2) + np.power(TrP[i-1,j+1],2)
			V2a = V2a + a
			V3a = V3a + a*TrP[i-1,(j-1):(j+1)]

		D = D + d1/r

		V2 = np.array([V2, V2a/r])
		b1 = b1 + V3a/r

	V1 = np.mean(TrP)
	
	p = np.size(V1,2)

	e1 = 0
	E = np.zeros(3)
	f1 = np.array([0, 0, 0])
	F = np.array([0, 0, 0])

	for k in range(1,3,p):
		e1 = np.transpose(V1[(k-1):(k+1)]) * V1[(k-1):(k+1)]
		E = E + e1
		f1 = V2[(k-2)/3 + 1] * V1[(k-1):(k+1)]
		F = F + f1

	# Equation 5 of Gamage and Lasenby
	A = 2 * (D - E)
	B = np.transpose(b1 - F)
	U, S, V = np.linalg.svd(A)

	Cm_in = V * np.linalg.inv(S) * np.transpose(U) * B
	Cm_old = Cm_in + np.transpose(np.array([1, 1, 1]))

	while distance(np.transpose(Cm_old), np.transpose(Cm_in)) > 0.0000001:
		Cm_old = Cm_in
		sigma2 = np.array([])

		for j in range(1,3,c):
			marker = TrP[:,(j-1):(j+1)]
			Ukp = marker - np.transpose(Cm_in * np.ones(1,r))

			# Computation of u^2
			u2 = 0
			app = np.array([])

			for i in range(1,r):
				u2 = u2 + Ukp[i-1,:] * np.transpose(Ukp[i-1,:])
			
			u2 = u2/r

			# Computation of sigma
			sigmaP = 0

			for i in range(1,r):
				sigmaP = sigmaP + np.power((app[i-1] - u2),2)
			
			sigmaP = sigmaP/(4 * u2 * r)
			sigma2 = np.array([np.transpose(sigma2), sigmaP])

			# Computation of deltaB
			deltaB = 0

			for j in range(1,3,c):
				deltaB = deltaB + np.transpose(V1[(j-1):(j+1)]) - Cm_in
			
			deltaB = 2 * sigma2 * deltaB

			Bcorr = B - deltaB # Corrected term B
			
			# Iterative estimation of the centre of rotation
			U, S, V = np.linalg.svd(A)

			Cm_in = V * np.linalg.inv(S) * np.transpose(U) * Bcorr
		
		Cm = Cm_in

	return Cm