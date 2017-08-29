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
import cv2
from scipy.ndimage import sobel, generic_gradient_magnitude, gaussian_gradient_magnitude, gaussian_filter, gaussian_laplace
from mpl_toolkits.mplot3d import Axes3D
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
    ref = dicom.read_file(files[100])

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
    orientations = []

    # loop through all the DICOM files
    for filename in files:
        # read the file
        ds = dicom.read_file(filename)

        #get pixel spacing and origin information
        origins.append(ds.ImagePositionPatient) #[0,0,0] coordinates in real 3D space (in mm)
        pixel_spacings.append([ds.PixelSpacing[0], ds.PixelSpacing[1], ds.SliceThickness]) #Space between pixels in x, y, z directions
        orientations.append(ds.ImageOrientationPatient)
        # store the raw image data
        dcm[:, :, files.index(filename)] = ds.pixel_array
    return dcm, origins, pixel_spacings, orientations

def preprocessing(dcm, origins, pixel_spacings, orientations):
    # Gradient Magnitude
    magnitude, direction = calculate_gradient_magnitude(dcm)

    # IGM Bins/Histogram
    # bins, tuples = create_bins(dcm, magnitude)

    # Refinement
    # patient_positions = get_patient_position(dcm, origins, pixel_spacings, orientations)
    # patient_positions = np.load('./data/saved/world_coordinates.npy')

    # bins, tuples = refine_igm_histogram(patient_positions, bins, tuples, False)

    # LH Bins/Histogram
    create_lh_histogram(dcm, magnitude)

    # Similarity Matrix

def calculate_gradient_magnitude(dcm):
    print "calculating gradient magnitude"
    gradient_magnitude = []
    gradient_direction = []

    z_length = dcm.shape[2]

    for z in range(z_length):
        slice = dcm[:,:,z]
        slice = gaussian_filter(slice, sigma=1)

        # Sobel derivative filters
        imx = np.zeros(slice.shape)
        sobel(slice,0,imx)

        imy = np.zeros(slice.shape)
        sobel(slice,1,imy)

        magnitude = np.hypot(imx,imy)
        direction = np.arctan(imy,imx)
        gradient_magnitude.append(magnitude)
        gradient_direction.append(direction)

    gradient_magnitude = np.dstack(gradient_magnitude)
    gradient_direction = np.dstack(gradient_direction)

    return gradient_magnitude, gradient_direction

def create_bins(dcm,magnitude):
    print "separating voxels into bins"
    tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
    bins = []

    # DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
    # Iterate through magnitude and dcm arrays and separate voxels into bins
    dcm_it = np.nditer(dcm, flags=['multi_index'])
    while not dcm_it.finished:
        dcm_idx = dcm_it.multi_index
        x = dcm_idx[0]
        y = dcm_idx[1]
        z = dcm_idx[2]
        intensity = dcm[x,y,z]
        gradient_magnitude = magnitude[x,y,z]
        # intensity = dcm[x,y]
        # gradient_magnitude = magnitude[x,y]
        couple = (intensity, gradient_magnitude)
        try:
            index = tuples.index(couple)
            bins[index].append(dcm_idx)
        except ValueError:
            bins.append([dcm_idx])
            tuples.append(couple)
        dcm_it.iternext()
    print "finished"
    # print "saving bins and tuples"
    # np.save("./data/saved/single_slice_bins.npy", bins)
    # np.save("./data/saved/single_slice_tuples", tuples)
    return bins, tuples

def create_igm_histogram(dcm, magnitude):
    print "constructing IGM histogram"
    print "\t WARNING: This is very computationally expensive and will take a significant amount of time."
    tuples = []
    size = dcm.shape[0] * dcm.shape[1] * dcm.shape[2]
    x = np.empty(size, dtype=np.uint16)
    y = np.empty(size, dtype=np.uint16)
    idx = 0
    for x_val in range(0, dcm.shape[0]):
        for y_val in range(0,dcm.shape[1]):
            for z_val in range(0, dcm.shape[2]):
                x[idx] = dcm[x_val,y_val, z_val]
                y[idx] = magnitude[x_val,y_val, z_val]
                couple = (x[idx], y[idx])
                tuples.append(couple)
                idx += 1
    print "\tcreating color array"
    alphas = np.empty(size, dtype=np.float16)
    #count each instance of a tuple
    idx = 0
    for couple in tuples:
        count = tuples.count(couple)
        alphas[idx] = count
        idx += 1
    print "\tnormalizing alpha values"
    alphas = alphas/alphas.max(axis=0)
    print alphas
    rgba_colors = np.zeros((size,4))
    rgba_colors[:,0] = 1.0
    rgba_colors[:,3] = alphas
    plt.xlabel('intensity')
    plt.ylabel('gradient magnitude')
    plt.scatter(x, y, color=rgba_colors)
    # plt.show()

def get_patient_position(dcm, origins, pixel_spacings, orientations):
    """
        Image Space --> Anatomical (Patient) Space is an affine transformation
        using the Image Orientation (Patient), Image Position (Patient), and
        Pixel Spacing properties from the DICOM header
    """
    print "getting 3d voxel coordinates"
    # Iterate through all slices and stack results
    world_coordinates = np.empty((dcm.shape[0], dcm.shape[1],dcm.shape[2], 3))
    z_length = dcm.shape[2]
    affine_matrix = np.zeros((4,4), dtype=np.float32)
    for z in range(z_length):
        slice = dcm[:,:,z]

        # Get all data needed to construct affine matrix
        image_position = origins[z]
        image_orientation = orientations[z]
        image_orientation_x = image_orientation[:3]
        image_orientation_y = image_orientation[3:]
        pixel_spacing = pixel_spacings[z]
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
        slice_it = np.nditer(slice, flags=['multi_index'])
        while not slice_it.finished:
            idx = slice_it.multi_index
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
            world_coordinates[x,y,z] = coordinates
            slice_it.iternext()
        print "done with slice " , str(z+1)
    print "saving coordinates"
    np.save('./data/saved/world_coordinates.npy', world_coordinates)
    return world_coordinates

def refine_igm_histogram(patient_positions,bins, tuples, show_histogram):
    print "refining IGM histogram"
    # Variables
    bin_count = 0
    refined_bins = []
    refined_tuples = []
    index_counter = 0

    THRESHOLD = 0.4

    for bin in bins:
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0
        bin_size = len(bin)
        # Mean position of voxels
        for indices in bin:
            x = indices[0]
            y = indices[0]
            # Get patient position and add to sum
            x_sum += patient_positions[x,y,0]
            y_sum += patient_positions[x,y,1]
            z_sum += patient_positions[x,y,2]
        x_avg = x_sum/bin_size
        y_avg = y_sum/bin_size
        z_avg = z_sum/bin_size
        mean_position = [x_avg,y_avg,z_avg]

        # Variance of voxel positions
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
        # Thresholding
        if variance <= THRESHOLD:
            refined_bins.append(bin)
            indexed_couple = tuples[index_counter]
            couple = (indexed_couple[0], indexed_couple[1])
            refined_tuples.append(couple)
            bin_count += 1
        index_counter += 1

    if show_histogram:
        print "\tsetting up new histogram"
        x = np.empty(len(refined_bins), dtype=np.uint16)
        y = np.empty(len(refined_bins), dtype=np.uint16)
        idx = 0
        for couple in refined_tuples:
            x[idx] = couple[0]
            y[idx] = couple[1]
            idx+=1
        print "\tcreating color array"
        alphas = np.empty(len(refined_bins), dtype=np.float16)
        idx = 0
        for couple in refined_tuples:
            count = len(refined_bins[idx])
            alphas[idx] = count
            idx += 1
        print "\tnormalizing alpha values"
        alphas = alphas/alphas.max(axis=0)
        rgba_colors = np.zeros((len(refined_bins),4))
        rgba_colors[:,0] = 1.0
        rgba_colors[:,3] = alphas
        plt.xlabel('intensity')
        plt.ylabel('gradient magnitude')
        print "\tplotting..."
        plt.scatter(x, y, color=rgba_colors)
        plt.show()
    return refined_bins, refined_tuples

def create_lh_histogram(dcm, magnitude):
    print "constructing LH histogram"
    # Determine if voxels lie on boundary or not
    # threshold = 140
    # voxels = np.zeros(dcm.shape, dtype=np.float32)
    # print np.max(magnitude)
    # x_idx = 0
    # for lst in magnitude:
    #     y_idx = 0
    #     for item in lst:
    #         if item >= threshold:
    #             voxels[x_idx,y_idx] = item
    #         y_idx += 1
    #     x_idx += 1
    second_derivative = gaussian_laplace(dcm, sigma=1)
    # f = function that gives back gradient magnitude for a given intensity, position pair


def f(x,t):
    print ""
    # return magnitude x at given pixel position t

# TAKEN FROM:
# https://stackoverflow.com/questions/31206443/numpy-second-derivative-of-a-ndimensional-array
def hessian(x):
    """
    Calculate the hessian matrix with finite differences
    Parameters:
       - x : ndarray
    Returns:
       an array of shape (x.dim, x.ndim) + x.shape
       where the array[i, j, ...] corresponds to the second derivative x_ij
    """
    x_grad = np.gradient(x)
    hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype)
    for k, grad_k in enumerate(x_grad):
        # iterate over dimensions
        # apply gradient again to every component of the first derivative.
        tmp_grad = np.gradient(grad_k)
        for l, grad_kl in enumerate(tmp_grad):
            hessian[k, l, :, :] = grad_kl
    return hessian

# TAKEN FROM: http://www.math-cs.gordon.edu/courses/mat342/python/diffeq.py
def rk2a( f, x0, t ):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        x = rk2a(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.

    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """

    n = len( t )
    x = numpy.array( [ x0 ] * n )
    for i in xrange( n - 1 ):
        h = t[i+1] - t[i]
        k1 = h * f( x[i], t[i] ) / 2.0
        x[i+1] = x[i] + h * f( x[i] + k1, t[i] + h / 2.0 )

    return x


def construct_similarity_matrix(dcm, magnitude, bins, tuples, patient_positions):
    boundary = [92,195]
    non_boundary = [72,48]
    igm_similarity = np.zeros(dcm.shape, dtype=np.float32)

    # print igm_similarity.shape

def affinity_propagation(bins, tuples, data):
    print "clustering data"
    K = 5 # Starter value, will probably change


# def create_transfer_functions():

# def render():

if __name__ == '__main__':
    path = "data/1/"
    data,origins, pixel_spacings, orientations = collect_data(path)
    preprocessing(data, origins, pixel_spacings, orientations)
