
import numpy as np

def segment_orientation_V1V3(V1,V3):

	e1 = np.zeros((len(V1),3))
	e2 = np.zeros((len(V1),3))
	e3 = np.zeros((len(V1),3))

	for i in range(1,len(V1)):
		e1[i-1,:] = V1[i-1,:]/np.sqrt(np.dot(V1[i-1,:],V1[i-1,:]))
		e3[i-1,:] = V3[i-1,:]/np.sqrt(np.dot(V3[i-1,:],V3[i-1,:]))

	for i in range(1,len(V1)):
		e2[i-1,:] = np.cross(e3[i-1,:],e1[i-1,:])

	for i in range(1,len(V1)):
		e3[i-1,:] = np.cross(e1[i-1,:],e2[i-1,:])

	return e1, e2, e3