#!/usr/bin/env python3
import os
import nibabel as nib
from nibabel.gifti import GiftiDataArray, GiftiImage
import numpy as np
import argparse
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation


def calc_connectivity_list(surf):
    """
    Precompute the connectivity list for the surface mesh.
    """
    n_vertices = surf['vertices'].shape[0]
    connectivity_list = [[] for _ in range(n_vertices)]
    
    for face in surf['faces']:
        connectivity_list[face[0]].append(face[1])
        connectivity_list[face[0]].append(face[2])
        connectivity_list[face[1]].append(face[0])
        connectivity_list[face[1]].append(face[2])
        connectivity_list[face[2]].append(face[0])
        connectivity_list[face[2]].append(face[1])
    
    return connectivity_list

def calc_vertex_area(surf):
    """
    Calculate the vertex area for a given surface mesh.
    
    Parameters
    ----------
    surf : dict
        Dictionary containing 'vertices' and 'faces' of the surface mesh.
    
    Returns
    -------
    area : np.ndarray
        Array of vertex areas.
    """
    vertices = surf['vertices']
    faces = surf['faces']
    
    # Calculate the area of each face using the cross product
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]
    
    face_areas = np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1) / 2.0
    
    # Initialize vertex area array
    area = np.zeros(vertices.shape[0], dtype=np.float32)
    
    # Accumulate face areas to corresponding vertices
    for i in range(faces.shape[0]):
        area[faces[i]] += face_areas[i] / 3.0  # Each vertex contributes equally to the face area
    
    return area

def find_sorted_vertex_clusters(sorted_vertices, connectivity_list, area = None, min_cluster_size=100):
    """
    Find connected local clusters of vertices based on their pre-sorted values and connectivity.

    This function iterates over the sorted vertices in descending order (effectively lowering 
    the threshold) and tracks clusters of connected vertices, merging clusters as they become connected.
    Before two clusters, each larger than min_cluster_size, merge, they are considered final clusters.
    This ensures each cluster represents the largest extent of locally increased values without absorbing
    neighboring clusters. The function returns the final clusters associated with the sorted vertices.

    Parameters
    ----------
    sorted_vertices : np.ndarray
        Array of vertex indices sorted by their values in descending order.
    connectivity_list : list of lists
        Precomputed connectivity list where each index corresponds to a vertex and contains a list of connected vertices.
    area : np.ndarray, optional
        Array of vertex areas. If provided, it can be used to determine the size of clusters based on area instead of vertex count.
        Default is None, in which case the size is determined by the number of vertices.
    min_cluster_size : int, optional
        Minimum size of a cluster to be considered final. Default is 100.

    Returns
    -------
    clusters : np.ndarray
        Array of cluster indices corresponding to the sorted vertices, where each index indicates the cluster to which the vertex belongs.
    """
    if area is None:
        # if no vertex areas are provided, use the number of vertices as the size
        area = np.ones(len(sorted_vertices), dtype=np.float32)

    # determine the largest vertex index, so that we can initialize the clusters array
    max_idx = np.max(sorted_vertices)
    # initialize array to grow clusters
    clusters_current = np.zeros(max_idx+1,dtype=np.uint32)
    # initialize array to store final clusters
    clusters_final = np.zeros(max_idx+1, dtype=np.uint32)

    # cluster numbering starts from 2, 1 is reserved for clusters that grew out of final clusters
    # 0 is reserved for unprocessed vertices
    next_cluster_idx = 2

    # iterate over vertices sorted by their values in descending order (effectively lowering the threshold)
    for vertex in sorted_vertices:  
        # before any potential merging, assign a new cluster index to the current vertex
        clusters_current[vertex] = next_cluster_idx
        next_cluster_idx += 1

        # find all connected vertices and check wheter they already belong to a cluster
        connected_vertices = np.minimum(np.array(connectivity_list[vertex]),max_idx)
        connected_cluster_idcs = {c for c in clusters_current[connected_vertices] if c != 0}
        
        # iteratively merge with all connected clusters
        cluster_idx_1 = clusters_current[vertex]
        for cluster_idx_2 in connected_cluster_idcs:
            # find all vertices belonging to the each cluster
            cluster_vertices_1 = np.where(clusters_current == cluster_idx_1)[0]
            cluster_vertices_2 = np.where(clusters_current == cluster_idx_2)[0]

            # determins current size of each cluster to be merged
            cluster_size_1 = np.sum(area[cluster_vertices_1])
            cluster_size_2 = np.sum(area[cluster_vertices_2])

            # if both clusters to be merged are larger than min_cluster_size, every not yet final cluster becomes a final cluster
            if cluster_size_1 >= min_cluster_size and cluster_size_2 >= min_cluster_size:
                if cluster_idx_1 != 1:
                    clusters_final[cluster_vertices_1] = cluster_idx_1
                if cluster_idx_2 != 1:
                    clusters_final[cluster_vertices_2] = cluster_idx_2
                merged_cluster_idx = 1
            # if one of the clusters is past final, the merged cluster also becomes past final
            elif cluster_idx_1 == 1 or cluster_idx_2 == 1:
                merged_cluster_idx = 1
            # otherwise get a new cluster index for the merged cluster
            else:
                merged_cluster_idx = next_cluster_idx
                next_cluster_idx += 1
            
            # merge the two clusters
            clusters_current[cluster_vertices_1] = merged_cluster_idx
            clusters_current[cluster_vertices_2] = merged_cluster_idx
        
            # update the cluster index of the first cluster to be merged
            cluster_idx_1 = merged_cluster_idx

    # return final clusters for all sorted vertices
    clusters = clusters_final[sorted_vertices]
    return clusters

def cluster_surface(metric, surf, area=None, 
                    min_cluster_size=100, min_abs_value=None, 
                    negative=False, output=None):
    """
    Cluster surface metric values into connected clusters of locally increased values.

    Parameters
    ----------
    metric : str or np.ndarray
        Path to the surface metric file (GIFTI format) or a numpy array of metric values.
    surf : str or dict
        Path to the surface mesh file (GIFTI format) or a dictionary with 'vertices' and 'faces'.
    area : str or np.ndarray, optional
        Path to the vertex area metric file (GIFTI format) or a numpy array of vertex areas.
        If not provided, area will be computed from the surface mesh.
    min_cluster_size : int, optional
        Minimum size of clusters to consider. Default is 100.
    min_abs_value : float, optional
        Minimum absolute value to consider for clustering. If not provided, uses the median of the absolute metric values.
    negative : bool, optional
        Whether to consider negative clusters. Default is False.
    output : str, optional
        Path to save the output clusters. If not provided, clusters will not be saved.

    Returns
    -------
    None
    """
    if isinstance(metric, str):
        metric_gii = nib.load(metric)
        metric = metric_gii.darrays[0].data
    
    if isinstance(surf, str):
        surf_gii = nib.load(surf)
        surf = dict(vertices=surf_gii.darrays[0].data,
                    faces=surf_gii.darrays[1].data)

    if area is None:
        area = calc_vertex_area(surf)
    elif isinstance(area, str):
        area_gii = nib.load(area)
        area = area_gii.darrays[0].data

    connectivity_list = calc_connectivity_list(surf)

    # first find all positive clusters using the minimum absolute value
    if min_abs_value is None:
        min_abs_value = np.nanmedian(np.abs(metric))

    # sort vertices by their metric values
    sorted_vertices = np.argsort(metric)[::-1]  # sort in descending order
    sorted_vertices_positive = sorted_vertices[metric[sorted_vertices] > min_abs_value]
    clusters_positive = find_sorted_vertex_clusters(sorted_vertices_positive, connectivity_list, area, min_cluster_size)

    if negative:
        # now find all negative clusters using the minimum absolute value
        sorted_vertices_negative = sorted_vertices[metric[sorted_vertices] < -min_abs_value][::-1]
        clusters_negative = find_sorted_vertex_clusters(sorted_vertices_negative, connectivity_list, area, min_cluster_size)

    # combine positive and optionally negative clusters as a surface metric, reindexing the clusters
    clusters = np.zeros_like(metric, dtype=np.uint32)
    unique_positive_clusters = np.unique(clusters_positive)
    for idx, cluster in enumerate(unique_positive_clusters[unique_positive_clusters > 0]):
        clusters[sorted_vertices_positive[clusters_positive == cluster]] = idx + 1
    if negative:
        unique_negative_clusters = np.unique(clusters_negative)
        for idx, cluster in enumerate(unique_negative_clusters[unique_negative_clusters > 0]):
            clusters[sorted_vertices_negative[clusters_negative == cluster]] = idx + 1 + len(unique_positive_clusters)

    # save the clusters to a GIFTI file if output is specified
    if output is not None:
        clusters_gii = GiftiImage()
        clusters_gii.add_gifti_data_array(GiftiDataArray(data=clusters, 
                                                              intent='NIFTI_INTENT_LABEL', 
                                                              datatype='NIFTI_TYPE_UINT32'))
        nib.save(clusters_gii, output)
    
    return clusters
    
def visualize_clusters(clusters, metric, surf):
    """
    Visualize the clustering result on the surface mesh. Works best with flattened surface as the 
    third dimension of the surface mesh is ignored for 2D visualization.
    
    Parameters
    ----------
    clusters : str or np.ndarray
        Path to the clusters file (GIFTI format) or a numpy array of cluster indices.
    metric : str or np.ndarray
        Path to the surface metric file (GIFTI format) or a numpy array of metric values.   
    surf : str or dict
        Path to the surface mesh file (GIFTI format) or a dictionary with 'vertices' and 'faces'.
    """
    if isinstance(metric, str):
        metric_gii = nib.load(metric)
        metric = metric_gii.darrays[0].data

    if isinstance(clusters, str):
        clusters_gii = nib.load(metric)
        clusters = clusters_gii.darrays[0].data

    if isinstance(surf, str):
        surf_gii = nib.load(surf)
        surf = dict(vertices=surf_gii.darrays[0].data,
                    faces=surf_gii.darrays[1].data)
        
    flatsurf = Triangulation(surf['vertices'][:, 0], surf['vertices'][:, 1], surf['faces'])

    # Mask zero values so they are not colored
    cmap = plt.cm.tab20([15,10,2,4,6,8,10,12,14,16,18,1,3,5,7,9,11,13,17,19])
    cmap = plt.cm.colors.ListedColormap(cmap)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].tripcolor(flatsurf, metric, shading='flat', cmap=plt.cm.hot)
    axes[0].axis('off')

    no_cluster = (clusters == 0)
    clusters = np.mod(clusters-1,19)+1
    clusters[no_cluster] = 0

    axes[1].tripcolor(flatsurf, clusters, shading='gouraud',vmin=0, vmax=20, cmap=cmap)
    axes[1].axis('off')

    # equal aspect ratio
    axes[0].set_aspect('equal')
    axes[1].set_aspect('equal')

    plt.tight_layout()
    plt.show()

    return fig, axes

def test_cluster_surface():
    """
    Test function to demonstrate the clustering on a sample surface and metric.
    """
    surf = os.path.join("/Users/dchaimow/data/test_cluster_surface",
                    "colin.cerebral.L.flat.164k_fs_LR.surf.gii")
    metric = os.path.join("/Users/dchaimow/data/test_cluster_surface",
                        "smoothed_glm_bold_fstat.func.gii")

    min_cluster_size=200    # minimum cluster size in mm^2 (default: 100 mm^2)
    min_abs_value=1         # minimum z-score (default: mean of absolute values)
    negative=False          # whether to consider negative clusters (default: False)

    clusters = cluster_surface(metric, surf,
                            min_cluster_size=min_cluster_size,
                            min_abs_value=min_abs_value,
                            negative=negative)
    visualize_clusters(clusters, metric, surf)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Cluster surface metric values into connected clusters of locally increased values.\n"
            "The algorithm iterates over vertices sorted by their metric values in descending order, "
            "tracking and merging clusters of connected vertices as they become adjacent. "
            "Clusters that reach the specified minimum size before merging are finalized, "
            "ensuring each cluster represents a locally maximal region. "
            "Only vertices with values above the minimum absolute value threshold are considered. "
            "Optionally, negative clusters can also be detected by processing vertices with values below the negative threshold. "
            "Final cluster assignments are saved as indices in a GIFTI metric file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("metric", type=str, help="Path to the surface metric file (GIFTI format).")
    parser.add_argument("surf", type=str, help="Path to the surface mesh file (GIFTI format).")
    parser.add_argument("output", type=str, help="Path to save the output clusters (GIFTI format).")
    parser.add_argument("--area", type=str, help="Path to the vertex area metric file (GIFTI format). If not provided, area will be calculated from the surface.")
    parser.add_argument("--min_cluster_size", type=int, default=100, help="Minimum size of clusters to consider.")
    parser.add_argument("--min_abs_value", type=float, default=None, help="Minimum absolute value to consider for clustering. If not provided, uses the median of absolute metric values.")
    parser.add_argument("--negative", action='store_true', help="Whether to consider negative clusters.")
    
    args = parser.parse_args()
    
    cluster_surface(args.metric, args.surf, area=args.area, 
                   min_cluster_size=args.min_cluster_size, 
                   min_abs_value=args.min_abs_value, 
                   negative=args.negative, output=args.output)

