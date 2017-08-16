import numpy as np
import dicom
import os
import matplotlib.pyplot as plt
from glob import glob
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import scipy.ndimage
from skimage import morphology
from skimage import measure
from skimage.transform import resize
from sklearn.cluster import KMeans
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.tools import FigureFactory as FF
from plotly.graph_objs import *
init_notebook_mode(connected=True)


def make_mesh(image, threshold=-300, step_size=1):

    print "Transposing surface"
    p = image.transpose(2,1,0)

    print "Calculating surface"
    verts, faces, norm, val = measure.marching_cubes(p, threshold, step_size=step_size, allow_degenerate=True)
    return verts, faces

def plotly_3d(verts, faces):
    x,y,z = zip(*verts)

    print "Drawing"

    # Make the colormap single color since the axes are positional not intensity.
#    colormap=['rgb(255,105,180)','rgb(255,255,51)','rgb(0,191,255)']
    colormap=['rgb(236, 236, 212)','rgb(236, 236, 212)']

    fig = FF.create_trisurf(x=x,
                        y=y,
                        z=z,
                        plot_edges=False,
                        colormap=colormap,
                        simplices=faces,
                        backgroundcolor='rgb(64, 64, 64)',
                        title="Interactive Visualization")
    iplot(fig)

def plt_3d(verts, faces):
    print "Drawing"
    x,y,z = zip(*verts)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces], linewidths=0.05, alpha=1)
    face_color = [1, 1, 0.9]
    mesh.set_facecolor(face_color)
    ax.add_collection3d(mesh)

    ax.set_xlim(0, max(x))
    ax.set_ylim(0, max(y))
    ax.set_zlim(0, max(z))
    ax.set_axis_bgcolor((0.7, 0.7, 0.7))
    plt.show()

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
v, f = make_mesh(dcm, 350)
plt_3d(v, f)
