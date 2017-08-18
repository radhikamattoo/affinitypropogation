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
    ref = dicom.read_file(files[127])
    return ref.pixel_array, [], [], ref


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
    # origins = []
    # pixel_spacings = []
    #
    # # loop through all the DICOM files
    # for filename in files:
    #     # read the file
    #     ds = dicom.read_file(filename)
    #
    #     #get pixel spacing and origin information
    #     origins.append(ds.ImagePositionPatient) #[0,0,0] coordinates in real 3D space (in mm)
    #     pixel_spacings.append([ds.PixelSpacing[0], ds.PixelSpacing[1], ds.SliceThickness]) #Space between pixels in x, y, z directions
    #
    #     # store the raw image data
    #     dcm[:, :, files.index(filename)] = ds.pixel_array
    # return dcm, origins, pixel_spacings

##############################################################################
#  Gradient Magnitude, Bins, Mean & Variance
##############################################################################
# Preprocessing includes constructing bins and calculating mean & variance
def preprocessing(dcm, origins, pixel_spacings, reference):
    print "calculating gradient magnitude"
    # Sobel filter for edge detection
    magnitude = generic_gradient_magnitude(dcm, sobel)

    # plt.subplot(1,2,1)
    # plt.imshow(magnitude)
    # plt.title('1st derivative')
    # plt.subplot(1,2,2)
    # plt.imshow(hess)
    # plt.title('2nd derivative')
    # plt.show()
    # sys.exit(0)

    # tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
    # bins = []

    # DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
    # Iterate through magnitude and dcm arrays and separate voxels into bins
    # dcm_it = np.nditer(dcm, flags=['multi_index'])
    # print "separating voxels into bins"
    # while not dcm_it.finished:
    #     dcm_idx = dcm_it.multi_index
    #     x = dcm_idx[0]
    #     y = dcm_idx[1]
    #     # z = dcm_idx[2]
    #     # intensity = dcm[x,y,z]
    #     # gradient_magnitude = magnitude[x,y,z]
    #     intensity = dcm[x,y]
    #     gradient_magnitude = magnitude[x,y]
    #     couple = (intensity, gradient_magnitude)
    #     try:
    #         index = tuples.index(couple)
    #         bins[index].append(dcm_idx)
    #     except ValueError:
    #         bins.append([dcm_idx])
    #         tuples.append(couple)
    #     dcm_it.iternext()
    # print "finished, saving bins and tuples"
    # np.save("./data/saved/single_slice_bins.npy", bins)
    # np.save("./data/saved/single_slice_tuples", tuples)
    # return bins, tuples
    print "loading bins and tuples"
    bins = np.load("./data/saved/single_slice_bins.npy")
    tuples = np.load("./data/saved/single_slice_tuples.npy")
    # print bins.shape
    # print tuples.shape
    # sys.exit(0)
    # x = np.empty(65536, dtype=np.uint16)
    # y = np.empty(65536, dtype=np.uint16)
    # print loaded_tuples.shape
    # for idx in range(0,65536):
    #     x[idx] = loaded_tuples[idx,0]
    #     y[idx] = loaded_tuples[idx,1]


    # print "attempting IGM histogram"
    # tuples = []
    # x = np.empty(65536, dtype=np.uint16)
    # y = np.empty(65536, dtype=np.uint16)
    # idx = 0
    # for x_val in range(0, 256):
    #     for y_val in range(0,256):
    #         x[idx] = dcm[x_val,y_val]
    #         y[idx] = magnitude[x_val,y_val]
    #         couple = (x[idx], y[idx])
    #         tuples.append(couple)
    #         idx += 1
    # print "creating color array"
    # alphas = np.empty(65536, dtype=np.float16)
    # #count each instance of a tuple
    # idx = 0
    # for couple in tuples:
    #     count = tuples.count(couple)
    #     alphas[idx] = count
    #     idx += 1
    # print "normalizing alpha values"
    # alphas = alphas/alphas.max(axis=0)
    # print alphas
    # rgba_colors = np.zeros((65536,4))
    # rgba_colors[:,0] = 1.0
    # rgba_colors[:,3] = alphas
    # plt.xlabel('intensity')
    # plt.ylabel('2nd derivative')
    # # plt.scatter(x, y, color=rgba_colors)
    # plt.scatter(x, y)
    # plt.show()


    print "constructing 3D positions of voxels"
    patient_positions = get_patient_position(dcm,origins,pixel_spacings,reference)

    print "refining IGM histogram"
    bin_count = 0
    refined_bins = []
    refined_tuples = []
    THRESHOLD = 0.4
    index_counter = 0
    for bin in bins:
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0
        bin_size = len(bin)
        #mean position of voxels
        for indices in bin:
            x = indices[0]
            y = indices[0]
            #get patient position
            #add to sum
            x_sum += patient_positions[x,y,0]
            y_sum += patient_positions[x,y,1]
            z_sum += patient_positions[x,y,2]
        x_avg = x_sum/bin_size
        y_avg = y_sum/bin_size
        z_avg = z_sum/bin_size
        mean_position = [x_avg,y_avg,z_avg]

        #variance of voxel positions
        variance = 0
        for indices in bin:
            x = indices[0]
            y = indices[0]
            position_x = patient_positions[x,y,0]
            position_y = patient_positions[x,y,1]
            position_z = patient_positions[x,y,2]

            difference_x = np.square(position_x - mean_position[0])
            difference_y = np.square(position_y - mean_position[1])
            difference_z = np.square(position_z - mean_position[2])

            variance += np.sqrt(difference_x + difference_y + difference_z)
        variance /= bin_size
        if variance <= THRESHOLD:
            refined_bins.append(bin)
            indexed_couple = tuples[index_counter]
            couple = (indexed_couple[0], indexed_couple[1])
            refined_tuples.append(couple)
            bin_count += 1
        index_counter += 1

    print len(refined_bins), len(refined_tuples)
    print "refined IGM histogram"
    x = np.empty(len(refined_bins), dtype=np.uint16)
    y = np.empty(len(refined_bins), dtype=np.uint16)
    idx = 0
    for couple in refined_tuples:
        x[idx] = couple[0]
        y[idx] = couple[1]
        idx+=1
    print "creating color array"
    alphas = np.empty(len(refined_bins), dtype=np.float16)
    idx = 0
    for couple in refined_tuples:
        count = len(refined_bins[idx])
        print count
        alphas[idx] = count
        idx += 1
    print "normalizing alpha values"
    alphas = alphas/alphas.max(axis=0)
    rgba_colors = np.zeros((len(refined_bins),4))
    rgba_colors[:,0] = 1.0
    rgba_colors[:,3] = alphas
    plt.xlabel('intensity')
    plt.ylabel('gradient magnitude')
    plt.scatter(x, y, color=rgba_colors)
    # plt.scatter(x, y)
    plt.show()
    sys.exit(0)




def get_patient_position(dcm, origins, pixel_spacings, dicom_object):
    """
        Image Space --> Anatomical (Patient) Space is an affine transformation
        using the Image Orientation (Patient), Image Position (Patient), and
        Pixel Spacing properties from the DICOM header
    """
    affine_matrix = np.zeros((4,4), dtype=np.float32)

    # Get all data needed to construct affine matrix
    image_position = dicom_object.ImagePositionPatient
    image_orientation = dicom_object.ImageOrientationPatient
    image_orientation_x = image_orientation[:3]
    image_orientation_y = image_orientation[3:]
    pixel_spacing = dicom_object.PixelSpacing
    pixel_spacing_x = pixel_spacing[0]
    pixel_spacing_y = pixel_spacing[1]

    # print image_position
    # print image_orientation_x, image_orientation_y
    # print pixel_spacing_x, pixel_spacing_y

    affine_matrix[0,0] = image_orientation_y[0] * pixel_spacing_x
    affine_matrix[0,1] = image_orientation_x[0] * pixel_spacing_y
    affine_matrix[0,3] = image_position[0]
    affine_matrix[1,0] = image_orientation_y[1] * pixel_spacing_x
    affine_matrix[1,1] = image_orientation_x[1] * pixel_spacing_y
    affine_matrix[1,3] = image_position[1]
    affine_matrix[2,0] = image_orientation_y[2] * pixel_spacing_x
    affine_matrix[2,1] = image_orientation_x[2] * pixel_spacing_y
    affine_matrix[2,3] = image_position[2]
    affine_matrix[3,3] = 1

    # Iterate through dcm pixels and perform affine transformation
    world_coordinates = np.empty((256,256,3))
    dcm_it = np.nditer(dcm, flags=['multi_index'])
    while not dcm_it.finished:
        idx = dcm_it.multi_index
        x = idx[0]
        y = idx[1]
        # Construct 4x1 vector
        vector = np.zeros((4,1), dtype=np.float32)
        vector[0,0] = x
        vector[1,0] = y
        vector[3,0] = 1
        # Perform affine transformation
        coordinates = np.matmul(affine_matrix, vector)
        coordinates = np.delete(coordinates, 3, axis=0)
        coordinates = np.transpose(coordinates)
        world_coordinates[x,y] = coordinates
        dcm_it.iternext()
    # print world_coordinates
    return world_coordinates

###################################################
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
    data, origins, pixel_spacings, reference = collect_data(path)
    preprocessing(data, origins, pixel_spacings, reference)
