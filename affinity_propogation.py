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
def collect_data(data_path):
    print "collecting data"
    files = []  # create an empty list
    for dirName, subdirList, fileList in os.walk(data_path):
        for filename in fileList:
            if ".dcm" in filename:
                files.append(os.path.join(dirName,filename))
    # Get reference file
    ref = dicom.read_file(files[0])
    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    pixel_dims = (int(ref.Rows), int(ref.Columns), len(files))

    # Load spacing values (in mm)
    pixel_space = (float(ref.PixelSpacing[0]), float(ref.PixelSpacing[1]), float(ref.SliceThickness))

    x = np.arange(0.0, (pixel_dims[0]+1)*pixel_space[0], pixel_space[0])
    y = np.arange(0.0, (pixel_dims[1]+1)*pixel_space[1], pixel_space[1])
    z = np.arange(0.0, (pixel_dims[2]+1)*pixel_space[2], pixel_space[2])

    # This will become the intensity values
    dcm = np.zeros(pixel_dims, dtype=ref.pixel_array.dtype)

    origins = []
    pixel_spacings = []

    # loop through all the DICOM files
    for filename in files:
        # read the file
        ds = dicom.read_file(filename)

        #get pixel spacing and origin information
        origins.append(ds.ImagePositionPatient) #[0,0,0] coordinates in real 3D space (in mm)
        pixel_spacings.append([ds.PixelSpacing[0], ds.PixelSpacing[1], ds.SliceThickness]) #Space between pixels in x, y, z directions

        # store the raw image data
        dcm[:, :, files.index(filename)] = ds.pixel_array
    return dcm, origins, pixel_spacings

##############################################################################
#  Gradient Magnitude, Bins, Mean & Variance
##############################################################################
# Preprocessing includes constructing bins and calculating mean & variance
def preprocessing(dcm, origins, pixel_spacings):
    print "calculating gradient magnitude"
    # Sobel filter for edge detection
    magnitude = generic_gradient_magnitude(data, sobel)
    # tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
    # bins = []
    #
    # # DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
    # # Iterate through magnitude and dcm arrays and separate voxels into bins
    # dcm_it = np.nditer(dcm, flags=['multi_index'])
    # print "separating voxels into bins"
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
    # print "finished, saving bins and tuples"
    # np.save("bins", bins)
    # np.save("tuples", tuples)
    # return bins, tuples
    print "loading bins and tuples"
    loaded_bins = np.load("./data/saved/bins.npy")
    loaded_tuples = np.load("./data/saved/tuples.npy")
    print "attempting IGM histogram"
    x = np.empty(loaded_tuples.shape[0], dtype=np.uint16)
    y = np.empty(loaded_tuples.shape[0], dtype=np.uint16)
    for idx, couple in enumerate(loaded_tuples):
        x[idx] = couple[0]
        y[idx] = couple[1]
    plt.xlabel('intensity')
    plt.ylabel('gradient magnitude')
    plt.scatter(x, y, alpha=0.09, s=1)
    plt.show()
    sys.exit(0)

    # print "constructing 3D positions of voxels"
    # positions = construct_3d(dcm, origins, pixel_spacings)

def construct_3d(voxels, origins, pixel_spacings):
    positions_3d = np.empty(voxels.shape, dtype=np.float32)


##############################################################################
# Refinement through removal of noisy 'bins'
##############################################################################
# def refinement(bins):

##############################################################################
# Affinity Propogation
##############################################################################
# def affinity_propogation(data):


##############################################################################
# Transfer Function Generation
##############################################################################
# def create_transfer_functions():

##############################################################################
# Volume Rendering
##############################################################################
# def render():
if __name__ == '__main__':
    path = "data/1/"
    data, origins, pixel_spacings = collect_data(path)
    preprocessing(data, origins, pixel_spacings)
