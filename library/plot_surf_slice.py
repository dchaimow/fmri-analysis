import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from voxeldepths_from_surfaces import load_fs_surf_in_grid


def plot_surf_slice(surf, z_slice, color='k'):
    """
    Plot a slice of a surface at z=z_slice
    surf: a dictionary with keys 'vertices' and 'tris'
    z_slice: the z value of the slice
    color: the color of the slice
    TODO: Implement slices for all three axes.
    """

    # Get vertices and faces from the surface
    vertices = surf['vertices']
    faces = surf['tris']

    # For all vertices, calculate on which side of the slice they are
    d = vertices[:,2] - z_slice
    out_vertices_mask = d > 0

    # For all faces, calculate how many vertices are on each side of the slice
    #  (i.e. how many vertices are in the mask)
    faces_out_vertices_count = np.sum(out_vertices_mask[faces] == True, axis=1)
    
    # If the number of outside vertices is 1 or 2, the face is cut by the slice
    cut_faces_mask = np.logical_or(faces_out_vertices_count == 1 ,
                                faces_out_vertices_count == 2)
    cut_faces = faces[cut_faces_mask]
    n_cut_faces = cut_faces.shape[0]
    
    # Flip the in/out assignemnt of the face vertices with 1 outside vertex
    #  because the math works the same way
    cut_faces_out_vertices = out_vertices_mask[cut_faces]
    to_flip = np.sum(cut_faces_out_vertices, axis=1) == 1
    cut_faces_out_vertices[to_flip] = np.logical_not(cut_faces_out_vertices[to_flip])

    # Sort cut face vertices into the one inside and the two outside
    single_side_indices = np.where(~cut_faces_out_vertices)[1]
    double_side_indices = np.where(cut_faces_out_vertices)
    first_double_side_indices = double_side_indices[1][::2]
    second_double_side_indices = double_side_indices[1][1::2]

    # Calculate the intersection points (c1,c2) of the slice with lines connecting
    #  the vertices of the cut faces (a to b1) and (a to b2)
    a = vertices[cut_faces[np.arange(n_cut_faces),single_side_indices], :]
    b1 = vertices[cut_faces[np.arange(n_cut_faces),first_double_side_indices], :]
    b2 = vertices[cut_faces[np.arange(n_cut_faces),second_double_side_indices], :]    
    t1 = (z_slice - a[:, 2]) / (b1[:, 2] - a[:, 2])
    t2 = (z_slice - a[:, 2]) / (b2[:, 2] - a[:, 2])
    c1 = a + (b1 - a) * t1[:, np.newaxis]
    c2 = a + (b2 - a) * t2[:, np.newaxis]

    plt.plot(
        np.stack([c1[:, 0], c2[:, 0]], axis=0),
        np.stack([c1[:, 1], c2[:, 1]], axis=0),
        f'{color}-', lw=0.5
    )


