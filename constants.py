#%%############################################################################
# Chose the parameters for the 1D readout geometry
###############################################################################
STRIP_WIDTH = 0.2  # (mm). The typical strip width is 0.2mm
PITCH = 100 / 256  # (mm). The typical strip pitch is 0.39mm (100mm/256strips)

###############################################################################
# Defining the strip certer positions
# The region R is where we are going to calculate the charge collected by the strips is:
#   first_strip_pos< R< last_strip_pos
# The strip_centers list define the position of the strip center
###############################################################################
FIRST_STRIP_POSITION = -5  # (mm)
LAST_STRIP_POSITION = 5  # (mm)

###############################################################################
# Defining positions of the electron cloud centers
# The region R2 is where we are going to calculate the charge collected by the strips is:
#   first_cloud_pos< R2 <last_cloud_pos
# The electron_cloud_centers list define the position where the electron cloud centers are going to be
###############################################################################
FIRST_CLOUD_POSITION = -2 * PITCH  # (mm)
LAST_CLOUD_POSITION = 2 * PITCH  # (mm)

###############################################################################
# Parameters for the clusterization
# the seed  value must be less than the maximum charge that a strip can collect considering the electron cloud choosen
# A charge value collected by a strip that is lower than the threshold is not considered part of an electron cloud by the algorithm
###############################################################################
SEED = 0.08
THRESHOLD = 0

STD_DEVIATION_OF_THE_NOISE = 0.01
NUMBER_OF_ELECTRON_CLOUDS = 1000
