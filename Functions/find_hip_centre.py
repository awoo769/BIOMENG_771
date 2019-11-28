
import numpy as np
from read_trc import read_trc
from segment_orientation_V1V3 import segment_orientation_V1V3
from calc_HJC import calc_HJC

mkr_data = read_trc()

LASI = mkr_data["Data"]["Markers"]["LAsis"]["All"]
RASI = mkr_data["Data"]["Markers"]["RAsis"]["All"]
LPSI = mkr_data["Data"]["Markers"]["LPsis"]["All"]
RPSI = mkr_data["Data"]["Markers"]["RPsis"]["All"]

TH1 = mkr_data["Data"]["Markers"]["RThighSuperior"]["All"]
TH2 = mkr_data["Data"]["Markers"]["RThighInferior"]["All"]
TH3 = mkr_data["Data"]["Markers"]["RThighLateral"]["All"]

# Define sacrum marker
SACR = (LPSI + RPSI)/2

# Define the pelvis origin
origin_pelvis = (LASI + RASI)/2

# Calculate unit vectors of pelvis segment relative to global
e1_pelvis, e2_pelvis, e3_pelvis = segment_orientation_V1V3(origin_pelvis - SACR, RASI - LASI)

TH1_pelvis = []
TH2_pelvis = []
TH3_pelvis = []

# Segment_orientation_V1V3 is good. Continue with writing findhipcentre to test the other 2.
# Transform the marker points into the pelvis coordinate system
for i in range(len(e1_pelvis)):
	rotation_matrix = np.array([e1_pelvis[i,:], e2_pelvis[i,:], e3_pelvis[i,:]])
	origin_pelvis_vector = np.transpose(origin_pelvis[i,:])

	TH1_pelvis.append(np.matmul(rotation_matrix,TH1[i,:]) - np.matmul(rotation_matrix,origin_pelvis_vector))
	TH2_pelvis.append(np.matmul(rotation_matrix,TH2[i,:]) - np.matmul(rotation_matrix,origin_pelvis_vector))
	TH3_pelvis.append(np.matmul(rotation_matrix,TH3[i,:]) - np.matmul(rotation_matrix,origin_pelvis_vector))

TH1_pelvis = np.array(TH1_pelvis)
TH2_pelvis = np.array(TH2_pelvis)
TH3_pelvis = np.array(TH3_pelvis)

thigh_markers_in_pelvis_cs = np.concatenate((TH1_pelvis, TH2_pelvis, TH3_pelvis), axis=1)

HJC_location = calc_HJC(thigh_markers_in_pelvis_cs)
