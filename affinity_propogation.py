#!/usr/bin/env python
"""
==================================================
Affinity Propogation for Improved Volume Rendering
==================================================
@author: Radhika Mattoo, http://www.github.com/radhikamattoo
@organization: MediVis, Inc.
@contact: rm3485@nyu.edu

Reference:
Tianjin Zhang et al., "A Clustering-Based Automatic Transfer Function Design for
Volume Visualization", Hindawi Sept. 2016
"""
print(__doc__)
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import medpy
import scipy
from scipy import ndimage
from medpy.io import load
from medpy.features.intensity import intensities

##############################################################################
# Data Collection & Preprocessing
##############################################################################
data_path = "data/CR-MONO1-10-chest"
image_data, image_header = load(data_path)
image_data = np.rot90(image_data,k=3)
print image_data.shape, image_data.dtype
intensity = intensities(image_data)
# print intensities(image_data)

# Calculate gradient magnitude and angle using Sobel filters
gx = ndimage.sobel(image_data, axis=0, mode='constant')
gy = ndimage.sobel(image_data, axis=1, mode='constant')
mag = np.sqrt(np.square(gx) + np.square(gy))
angle = np.arctan(gy/gx)

# Plot gradient
# plt.close("all")
# plt.figure()
# plt.suptitle("Image, and it gradient along each axis")
# ax = plt.subplot("131")
# ax.axis("off")
# ax.imshow(mag)
# ax.set_title("magnitude")
#
# ax = plt.subplot("132")
# ax.axis("off")
# ax.imshow(gx)
# ax.set_title("gx")
#
# ax = plt.subplot("133")
# ax.axis("off")
# ax.imshow(gy)
# ax.set_title("gy")
# plt.show()

# Construct IGM
plt.scatter(intensity, mag, alpha=0.006)
plt.show()

##############################################################################
# Variance
##############################################################################


##############################################################################
# Refinement through removal of noisy 'bins'
##############################################################################


##############################################################################
# Affinity Propogation
##############################################################################



##############################################################################
# Transfer Function Generation
##############################################################################



##############################################################################
# Volume Rendering
##############################################################################
