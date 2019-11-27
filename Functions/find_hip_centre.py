
from read_trc import read_trc
from segment_orientation_V1V3 import segment_orientation_V1V3

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

# Segment_orientation_V1V3 is good. Continue with writing findhipcentre to test the other 2.
# Transform the marker points into the pelvis coordinate system
#or i in range(1, len(e1_pelvis)):
#	rotation_matrix = np.array([e1_pelvis[i-1,:]])