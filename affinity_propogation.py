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
import nibabel as nib
import dicom
from scipy.ndimage import sobel, generic_gradient_magnitude
from medpy.io import load
from medpy.features.intensity import intensities, gaussian_gradient_magnitude

##############################################################################
# Data Collection & Preprocessing
##############################################################################
data_path = "data/OBJ_0001/"

files = []  # create an empty list
for dirName, subdirList, fileList in os.walk(data_path):
    for filename in fileList:
        if "IM_" in filename:
            files.append(os.path.join(dirName,filename))
# Get ref file
ref = dicom.read_file(files[0])
# Load dimensions based on the number of rows, columns, and slices (along the Z axis)
pixel_dims = (int(ref.Rows), int(ref.Columns), len(files))

# Load spacing values (in mm)
pixel_space = (float(ref.PixelSpacing[0]), float(ref.PixelSpacing[1]), float(ref.SliceThickness))

x = np.arange(0.0, (pixel_dims[0]+1)*pixel_space[0], pixel_space[0])
y = np.arange(0.0, (pixel_dims[1]+1)*pixel_space[1], pixel_space[1])
z = np.arange(0.0, (pixel_dims[2]+1)*pixel_space[2], pixel_space[2])

# This will become the intensity values used later on
dcm = np.zeros(pixel_dims, dtype=ref.pixel_array.dtype)

# loop through all the DICOM files
for filename in files:
    # read the file
    ds = dicom.read_file(filename)
    # store the raw image data
    dcm[:, :, files.index(filename)] = ds.pixel_array

# plt.figure(dpi=300)
# # plt.axes().set_aspect('equal', 'datalim')
# plt.set_cmap(plt.gray())
# plt.pcolormesh(x, y,z)
# plt.show()

# image_data = np.rot90(image_data,k=3)
# intensity = intensities(image_data)
# print image_data.shape, image_data.dtype
# print intensities(image_data)

##############################################################################
#  Gradient Magnitude & Bins
##############################################################################
# Sobel filters
print "calculating gradient magnitude"
magnitude = generic_gradient_magnitude(dcm, sobel)
tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
bins = []

# DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
# Iterate through magnitude and dcm arrays and separate voxels into bins
dcm_it = np.nditer(dcm, flags=['multi_index'])
print "separating voxels into bins"
count = 0
while not dcm_it.finished:
    dcm_idx = dcm_it.multi_index
    x = dcm_idx[0]
    y = dcm_idx[1]
    z = dcm_idx[2]
    intensity = dcm[x,y,z]
    gradient_magnitude = magnitude[x,y,z]
    couple = (intensity, gradient_magnitude)
    if couple not in tuples: #add a new bin
        # print "creating new bin"
        tuples.append(couple)
        bins.append({'indices':[dcm_idx], 'gradient_magnitude': gradient_magnitude, 'intensity': intensity})
    else: #we found the tuple, so add the index to an existing bin
        # print "adding to existing bin"
        for dictionary in bins:
            if dictionary['gradient_magnitude'] == gradient_magnitude and dictionary['intensity'] == intensity:
                dictionary['indices'].append(dcm_idx)
    count += 1
    if count % 1000000 == 0:
        print count
    dcm_it.iternext()
print "finished"
print len(bins)


##############################################################################
# Variance & Mean of Position
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
