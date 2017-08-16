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
from mpl_toolkits.mplot3d import Axes3D
import medpy
import scipy
import nibabel as nib
import dicom
from scipy.ndimage import sobel, generic_gradient_magnitude, gaussian_gradient_magnitude
from medpy.io import load

##############################################################################
# Data Collection & Preprocessing
##############################################################################
# print "collecting data"
# data_path = "data/1/"
# files = []  # create an empty list
# for dirName, subdirList, fileList in os.walk(data_path):
#     for filename in fileList:
#         if ".dcm" in filename:
#             files.append(os.path.join(dirName,filename))
# # Get reference file
# ref = dicom.read_file(files[0])
#
# # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
# pixel_dims = (int(ref.Rows), int(ref.Columns), len(files))
#
# # Load spacing values (in mm)
# pixel_space = (float(ref.PixelSpacing[0]), float(ref.PixelSpacing[1]), float(ref.SliceThickness))
#
# x = np.arange(0.0, (pixel_dims[0]+1)*pixel_space[0], pixel_space[0])
# y = np.arange(0.0, (pixel_dims[1]+1)*pixel_space[1], pixel_space[1])
# z = np.arange(0.0, (pixel_dims[2]+1)*pixel_space[2], pixel_space[2])
#
# # This will become the intensity values
# dcm = np.zeros(pixel_dims, dtype=ref.pixel_array.dtype)
#
# # loop through all the DICOM files
# for filename in files:
#     # read the file
#     ds = dicom.read_file(filename)
#     # store the raw image data
#     dcm[:, :, files.index(filename)] = ds.pixel_array
#
# ##############################################################################
# #  Gradient Magnitude & Bins
# ##############################################################################
# # Sobel filter for edge detection
# print "calculating gradient magnitude"
# # magnitude = generic_gradient_magnitude(dcm, sobel)
# magnitude = generic_gradient_magnitude(ref.pixel_array, sobel)
#
# tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
# bins = []

# DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
# Iterate through magnitude and dcm arrays and separate voxels into bins
# dcm_it = np.nditer(dcm, flags=['multi_index'])
# print "separating voxels into bins"
# count = 0
# while not dcm_it.finished:
#     dcm_idx = dcm_it.multi_index
#     x = dcm_idx[0]
#     y = dcm_idx[1]
#     z = dcm_idx[2]
#     intensity = dcm[x,y,z]
#     gradient_magnitude = magnitude[x,y,z]
#     couple = (intensity, gradient_magnitude)
#     if couple not in tuples: #add a new bin
#         # print "creating new bin"
#         tuples.append(couple)
#         bins.append({'indices':[dcm_idx], 'gradient_magnitude': gradient_magnitude, 'intensity': intensity})
#     else: #we found the tuple, so add the index to an existing bin
#         # print "adding to existing bin"
#         for dictionary in bins:
#             if dictionary['gradient_magnitude'] == gradient_magnitude and dictionary['intensity'] == intensity:
#                 dictionary['indices'].append(dcm_idx)
#     count += 1
#     if count % 100000 == 0:
#         print str(count) + "/" + str(13631488)
#     dcm_it.iternext()

# dcm_it = np.nditer(dcm, flags=['multi_index'])
# print "separating voxels into bins"
# count = 0
# while not dcm_it.finished:
#     dcm_idx = dcm_it.multi_index
#     x = dcm_idx[0]
#     y = dcm_idx[1]
#     z = dcm_idx[2]
#     intensity = dcm[x,y,z]
#     gradient_magnitude = magnitude[x,y,z]
#     couple = (intensity, gradient_magnitude)
#     try:
#         index = tuples.index(couple)
#         bins[index].append(dcm_idx)
#     except ValueError:
#         bins.append([dcm_idx])
#         tuples.append(couple)
#     dcm_it.iternext()
# print "finished, saving bins and tuples arrays"
# np.save("bins_10", bins)
# np.save("tuples_10", tuples)

##############################################################################
# Variance & Mean of Position
##############################################################################
print "finding variance"
loaded_bins = np.load("./bins_10.npy")
loaded_tuples = np.load("./tuples_10.npy")

x = np.empty((35092), dtype=np.uint16)
y = np.empty((35092), dtype=np.uint16)

for idx,couple in enumerate(loaded_tuples):
    intensity = couple[0]
    gradient_magnitude = couple[1]
    x[idx] = intensity
    y[idx] = gradient_magnitude
print np.max(x), np.max(y)
axes = plt.gca()
axes.set_xlim([0,350])
axes.set_ylim([0,260])
plt.scatter(x, y,  cmap='Binary')
plt.show()





# for bin in loaded_bins:
#     n = len(bin)
#     x_avg = 0
#     y_avg = 0
#     z_avg = 0
#     for voxel_idx in bin:
#         x_avg += voxel_idx[0]
#         y_avg += voxel_idx[1]
#         z_avg += voxel_idx[2]
#     x_avg /= n
#     y_avg /= n
#     z_avg /= n
#     avg = [x_avg,y_avg,z_avg]
#
#     x_variance = 0
#     y_variance = 0
#     z_variance = 0
#     for voxel_idx in bin:
#         x_variance += np.abs(voxel_idx[0] - avg[0])
#         y_variance += np.abs(voxel_idx[1] - avg[1])
#         z_variance += np.abs(voxel_idx[2] - avg[2])
#     x_variance /= n
#     y_variance /= n
#     z_variance /= n
#
#     variance = [x_variance, y_variance, z_variance]
#     print variance
#     sys.exit(0)
#

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
