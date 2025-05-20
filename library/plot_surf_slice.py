import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from voxeldepths_from_surfaces import load_fs_surf_in_grid


def plot_surf_slice(surf, z_slice, color='k'):
    # load surface file
    vertices = surf['vertices']
    faces = surf['tris']

    # 1. for all vertices, calculate on which side of the slice they are
    d = vertices[:,2] - z_slice
    out_vertices_mask = d > 0
    # 2. for all faces, calculate how many vertices are on each side of the slice
    #    (i.e. how many vertices are in the mask)
    faces_out_vertices_count = np.sum(out_vertices_mask[faces] == True, axis=1)
    # 3. If the number of outside vertices is 1 or 2, the face is cut by the slice
    cut_faces_mask = np.logical_or(faces_out_vertices_count == 1 ,
                                faces_out_vertices_count == 2)
    cut_faces = faces[cut_faces_mask]
    n_cut_faces = cut_faces.shape[0]
    # flip the in/out assignemnt of the face vertices with 1 outside vertex
    # because the math works the same way
    cut_faces_out_vertices = out_vertices_mask[cut_faces]
    to_flip = np.sum(cut_faces_out_vertices, axis=1) == 1
    cut_faces_out_vertices[to_flip] = np.logical_not(cut_faces_out_vertices[to_flip])

    # sort cat face vertices into the one inside and the two outside
    single_side_indices = np.where(~cut_faces_out_vertices)[1]
    double_side_indices = np.where(cut_faces_out_vertices)
    first_double_side_indices = double_side_indices[1][::2]
    second_double_side_indices = double_side_indices[1][1::2]

    a = vertices[cut_faces[np.arange(n_cut_faces),single_side_indices], :]
    b1 = vertices[cut_faces[np.arange(n_cut_faces),first_double_side_indices], :]
    b2 = vertices[cut_faces[np.arange(n_cut_faces),second_double_side_indices], :]
    c1 = np.zeros((n_cut_faces, 2))
    c2 = np.zeros((n_cut_faces, 2))
    c1[:, 0] = a[:, 0] + (b1[:, 0] - a[:, 0]) * (z_slice - a[:, 2]) / (b1[:, 2] - a[:, 2])
    c1[:, 1] = a[:, 1] + (b1[:, 1] - a[:, 1]) * (z_slice - a[:, 2]) / (b1[:, 2] - a[:, 2])

    c2[:, 0] = a[:, 0] + (b2[:, 0] - a[:, 0]) * (z_slice - a[:, 2]) / (b2[:, 2] - a[:, 2])
    c2[:, 1] = a[:, 1] + (b2[:, 1] - a[:, 1]) * (z_slice - a[:, 2]) / (b2[:, 2] - a[:, 2])

    plt.plot(
        np.stack([c1[:, 0], c2[:, 0]], axis=0),
        np.stack([c1[:, 1], c2[:, 1]], axis=0),
        f'{color}-', lw=0.5
    )
    # set the axis limits to the min and max of the vertices
    # set the aspect ratio to be equal
    plt.gca().set_aspect('equal', adjustable='box')


