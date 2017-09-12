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
import time
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy
import dicom
from scipy.ndimage import sobel, generic_gradient_magnitude, gaussian_gradient_magnitude, gaussian_filter, laplace, gaussian_laplace
from sklearn.cluster import AffinityPropagation

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
    pixel_spacings = (float(ref.PixelSpacing[0]), float(ref.PixelSpacing[1]), float(ref.SliceThickness))

    x = np.arange(0.0, (pixel_dims[0]+1)*pixel_spacings[0], pixel_spacings[0])
    y = np.arange(0.0, (pixel_dims[1]+1)*pixel_spacings[1], pixel_spacings[1])
    z = np.arange(0.0, (pixel_dims[2]+1)*pixel_spacings[2], pixel_spacings[2])

    # Row and column directional cosines
    orientation = ref.ImageOrientationPatient

    # This will become the intensity values
    dcm = np.zeros(pixel_dims, dtype=ref.pixel_array.dtype)

    origins = []

    # loop through all the DICOM files
    for filename in files:
        # read the file
        ds = dicom.read_file(filename)
        #get pixel spacing and origin information
        origins.append(ds.ImagePositionPatient) #[0,0,0] coordinates in real 3D space (in mm)

        # store the raw image data
        dcm[:, :, files.index(filename)] = ds.pixel_array
    return dcm, origins, pixel_spacings, orientation

def preprocessing(dcm, origins, pixel_spacing, orientation):
    # Gradient Magnitude
    magnitude, azimuthal, elevation = calculate_gradient_magnitude(dcm)

    # IGM Bins/Histogram
    # bins, tuples = create_bins(dcm, magnitude)
    # print "loading bins and tuples"
    bins = np.load("./data/saved/refined_3d_bins.npy")
    tuples = np.load("./data/saved/refined_3d_tuples.npy")

    # create_igm_histogram(dcm, magnitude)

    # Refinement
    # patient_positions = get_patient_position(dcm, origins, pixel_spacing, orientation)
    patient_positions = np.load('./data/saved/world_coordinates.npy')

    # bins, tuples = refine_igm_histogram(patient_positions, bins, tuples, True)

    # LH Bins/Histogram
    # create_lh_histogram(patient_positions,dcm, magnitude, azimuthal, elevation)

    # Similarity Matrix
    construct_similarity_matrix(dcm, magnitude, bins, tuples, patient_positions)

def calculate_gradient_magnitude(dcm):
    print "calculating gradient magnitude"
    gradient_magnitude = []
    gradient_direction = []

    gradx = np.zeros(dcm.shape)
    sobel(dcm,0,gradx)
    grady = np.zeros(dcm.shape)
    sobel(dcm,1,grady)
    gradz = np.zeros(dcm.shape)
    sobel(dcm,2,gradz)

    gradient = np.sqrt(gradx**2 + grady**2 + gradz**2)

    azimuthal = np.arctan2(grady, gradx)
    elevation = np.arctan(gradz,gradx)

    azimuthal = np.degrees(azimuthal)
    elevation = np.degrees(elevation)

    return gradient, azimuthal, elevation

    # z_length = dcm.shape[2]
    #
    # for z in range(z_length):
    #     slice = dcm[:,:,z]
    #     slice = gaussian_filter(slice, sigma=1)
    #
    #     # Sobel derivative filters
    #     imx = np.zeros(slice.shape)
    #     sobel(slice,0,imx)
    #
    #     imy = np.zeros(slice.shape)
    #     sobel(slice,1,imy)
    #
    #     magnitude = np.hypot(imx,imy)
    #     direction = np.arctan(imy,imx)
    #     gradient_magnitude.append(magnitude)
    #     gradient_direction.append(direction)
    #
    # gradient_magnitude = np.dstack(gradient_magnitude)
    # gradient_direction = np.dstack(gradient_direction)
    #
    # return gradient_magnitude, gradient_direction

def create_bins(dcm,magnitude):
    print "separating voxels into bins"
    print "\t WARNING: This is very computationally expensive and will take a significant amount of time."
    # TODO: make parallel?

    tuples = [] #will hold tuples of (intensity, magnitude) for bin separation
    bins = []

    # Keeping track of how long this takes
    # DICOM IMAGE INTENSITIES STORED IN 12 BITS, 0-4095 VALUES
    # Iterate through magnitude and dcm arrays and separate voxels into bins
    dcm_it = np.nditer(dcm, flags=['multi_index'])
    count = 0
    tic = time.clock() # Measure time it takes to run through 1/8 of voxels (around 8 million total)
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
        # Print time information
        if count % 100000 == 0:
            toc = time.clock()
            time_taken = toc - tic
            tic = time.clock() #reset timer
            print str(count/100000) + " took " + str(time_taken) + " seconds"
        count += 1
        dcm_it.iternext()
    print "saving bins and tuples"
    np.save("./data/saved/3d_bins.npy", bins)
    np.save("./data/saved/3d_tuples.npy", tuples)
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

def get_patient_position(dcm, origins, pixel_spacing, orientation):
    """
        Image Space --> Anatomical (Patient) Space is an affine transformation
        using the Image Orientation (Patient), Image Position (Patient), and
        Pixel Spacing properties from the DICOM header
    """
    print "getting 3d coordinates"

    world_coordinates = np.empty((dcm.shape[0], dcm.shape[1],dcm.shape[2], 3))
    affine_matrix = np.zeros((4,4), dtype=np.float32)

    rows = dcm.shape[0]
    cols = dcm.shape[1]
    num_slices = dcm.shape[2]

    image_orientation_x = np.array([ orientation[0], orientation[1], orientation[2]  ]).reshape(3,1)
    image_orientation_y = np.array([ orientation[3], orientation[4], orientation[5]  ]).reshape(3,1)
    pixel_spacing_x = pixel_spacing[0]

    # Construct affine matrix
    # Method from:
    # http://nipy.org/nibabel/dicom/dicom_orientation.html
    T_1 = origins[0]
    T_n = origins[num_slices-1]


    affine_matrix[0,0] = image_orientation_y[0] * pixel_spacing[0]
    affine_matrix[0,1] = image_orientation_x[0] * pixel_spacing[1]
    affine_matrix[0,3] = T_1[0]
    affine_matrix[1,0] = image_orientation_y[1] * pixel_spacing[0]
    affine_matrix[1,1] = image_orientation_x[1] * pixel_spacing[1]
    affine_matrix[1,3] = T_1[1]
    affine_matrix[2,0] = image_orientation_y[2] * pixel_spacing[0]
    affine_matrix[2,1] = image_orientation_x[2] * pixel_spacing[1]
    affine_matrix[2,3] = T_1[2]
    affine_matrix[3,3] = 1

    k1 = (T_1[0] - T_n[0])/ (1 - num_slices)
    k2 = (T_1[1] - T_n[1])/ (1 - num_slices)
    k3 = (T_1[2] - T_n[2])/ (1 - num_slices)

    affine_matrix[:3, 2] = np.array([k1,k2,k3])

    for z in range(num_slices):
        for r in range(rows):
            for c in range(cols):
                vector = np.array([r, c, 0, 1]).reshape((4,1))
                result = np.matmul(affine_matrix, vector)
                result = np.delete(result, 3, axis=0)
                result = np.transpose(result)
                world_coordinates[r,c,z] = result
        # print "Finished slice ", str(z)
    # np.save('./data/saved/world_coordinates_3d.npy', str(world_coordinates))
    return world_coordinates

def refine_igm_histogram(patient_positions,bins, tuples, show_histogram=False):
    print "refining IGM histogram"
    # Variables
    bin_count = 0
    refined_bins = []
    refined_tuples = []
    index_counter = 0

    THRESHOLD = 0.9

    for bin in bins:
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0
        bin_size = len(bin)
        # Mean position of voxels
        for indices in bin:
            x = indices[0]
            y = indices[1]
            z = indices[2]
            # Get patient position and add to sum
            patient_position = patient_positions[x,y,z]
            x_sum += patient_position[0]
            y_sum += patient_position[1]
            z_sum += patient_position[2]
        x_avg = x_sum/bin_size
        y_avg = y_sum/bin_size
        z_avg = z_sum/bin_size

        mean_position = [x_avg,y_avg,z_avg]

        # Variance of voxel positions
        variance = 0
        for indices in bin:
            x = indices[0]
            y = indices[1]
            z = indices[2]
            patient_position = patient_positions[x,y,z]
            position_x = patient_position[0]
            position_y = patient_position[1]
            position_z = patient_position[2]

            difference_x = np.square(position_x - mean_position[0])
            difference_y = np.square(position_y - mean_position[1])
            difference_z = np.square(position_z - mean_position[2])
            variance += np.sqrt(difference_x + difference_y + difference_z)
        variance /= bin_size
        # Thresholding
        if variance <= THRESHOLD:
            print variance
            refined_bins.append(bin)
            indexed_couple = tuples[index_counter]
            couple = (indexed_couple[0], indexed_couple[1])
            refined_tuples.append(couple)
            bin_count += 1
        index_counter += 1
    print "saving refined bins/tuples"
    np.save("./data/saved/refined_3d_bins_09.npy", refined_bins)
    np.save("./data/saved/refined_3d_tuples_09.npy", refined_tuples)

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

def create_lh_histogram(patient_positions, dcm, magnitude, azimuthal, elevation):
    print "constructing LH histogram"
    # Get 2nd derivative
    second_derivative = gaussian_filter(magnitude, sigma=1, order=1)

    # Determine if voxels lie on boundary or not (thresholding)
    dcm_it = np.nditer(dcm, flags=['multi_index'])
    THRESHOLD = 100
    REPLACER = -100
    count = 0
    while not dcm_it.finished:
        dcm_idx = dcm_it.multi_index
        x = dcm_idx[0]
        y = dcm_idx[1]
        z = dcm_idx[2]
        intensity = dcm[x,y,z]
        if intensity <= THRESHOLD:
            dcm[x,y,z] = REPLACER
            count += 1
    print count
    #Iterate through all thresholded voxels and integrate gradient field in
    # both directions using 2nd-order Runge-Kutta
    vox_it = voxels.nditer(voxels, flags=['multi_index'])
    while not vox_it.finished:
        # ???
        print "meh"

# Returns intensity value at intermediary voxel position
def trilinear_interpolation(dcm, vx,vy,vz):
    print "interpolating"

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
    shape =(bins.shape[0], bins.shape[0])
    igm_similarity = np.zeros(shape, dtype=np.float32)
    MIN = -1
    MAX = 100000

    outer_tuple_counter = 0

    # (intensity, gradient magnitude)
    for bin in bins:
        outer_intensity = tuples[outer_tuple_counter][0]
        outer_gradient_magnitude = tuples[outer_tuple_counter][1]
        inner_tuple_counter = 0
        for bin in bins:
            inner_intensity = tuples[inner_tuple_counter][0]
            inner_gradient_magnitude = tuples[inner_tuple_counter][1]

            # IGM Euclidean distance
            abs_intensity = np.square( np.absolute( outer_intensity - inner_intensity ) )
            abs_gm = np.square( np.absolute( outer_gradient_magnitude - inner_gradient_magnitude ) )
            euclidean = np.sqrt(abs_intensity + abs_gm)

            igm_similarity[outer_tuple_counter, inner_tuple_counter] = euclidean
            print str(euclidean) + " for bin index (" + str(outer_tuple_counter) + " , " + str(inner_tuple_counter) + ")"
            inner_tuple_counter += 1
        outer_tuple_counter += 1
    np.save('./data/saved/similarity_matrix.npy', igm_similarity)

def affinity_propagation(bins, tuples, data):
    print "clustering data"
    K = 5 # Starter value, will probably change



if __name__ == '__main__':
    path = "data/1/"
    data,origins, pixel_spacing, orientation = collect_data(path)
    preprocessing(data, origins, pixel_spacing, orientation)
