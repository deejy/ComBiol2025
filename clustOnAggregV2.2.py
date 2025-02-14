import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN
import argparse
import matplotlib.colors as mcolors
import csv

# Constants
DEFAULT_EPS = 54.0
DEFAULT_MIN_SAMPLES = 2
DEFAULT_SEED = 427
DEFAULT_DIAM = 54

def load_data(input_file):
    """Load particle coordinates from a text file."""
    try:
        if input_file.endswith('.conf'):
            # Read the .conf file, skipping the first two lines and using the second and third columns as x and y
            data = pd.read_csv(input_file, delim_whitespace=True, skiprows=2, usecols=[1, 2], header=None)
        else:
            # Load normally for other file types (e.g., .txt)
            data = pd.read_csv(input_file, delim_whitespace=True, usecols=[0, 1], header=None)
        data.columns = ['x', 'y']  # Rename columns to 'x' and 'y'
        
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        raise
        
def calculate_distances(coordinates):
    """Calculate pairwise distances between particles."""
    return squareform(pdist(coordinates))

def perform_clustering(distances, threshold, oligo):
    """Perform DBSCAN clustering on the distance matrix."""
    db = DBSCAN(eps=threshold, min_samples=oligo, metric='precomputed', algorithm='brute') #brute or auto is the same for "precomputed" distance matrix
    return db.fit_predict(distances)


def generate_colors_by_size(cluster_sizes):
    """Generate colors for each cluster size with similar tints but high contrast."""
    unique_sizes = sorted(cluster_sizes.unique())
    colors = []    
    
    # Define a base color; you can choose a base hue and adjust its brightness
    base_hues = np.linspace(0, 0.4, 5, endpoint=False)  # Generate unique hues
    base_hues2 = np.linspace(0.4, 0.5, 10, endpoint=False)  # Generate unique hues
    base_hues3 = np.linspace(0.5, 0.6, 400, endpoint=False)  # Generate unique hues    
    fused_hues = np.concatenate((base_hues, base_hues2, base_hues3 ))
    
#    base_hue = 0.5  # For example, a green color in the HSV color space
    for size in unique_sizes:
        # Generate a color based on the size
        color = mcolors.hsv_to_rgb((fused_hues[size-2], 1, 1))  # Convert HSV to RGB
        colors.append((size, color))
    return dict(colors)
            

def calculate_point_size(coordinates, width):
    """Calculate the size of the points based a fake plot."""
    
    # A fake plot to extract figure area
    plt.figure(figsize=(10, 10))
    plt.scatter(coordinates['x'],
                        coordinates['y'],
                        color="red", alpha=0.7, s=1)                      
    fig = plt.gcf()
    fig_width, fig_height = fig.get_size_inches()
    plt.clf() #Clear figure
    plt.close()
      
    # Calculate a base size for points based on the figure area
    fig_area = fig_width * fig_height  # Area in square inches
    # The scaling factor is adhoc and corresponds to a protein of 53 A diameter 
    ps_ref = fig_area * 8.5 
    #point size are homogeneous to surfaces
    rpoint_ref = np.sqrt(ps_ref/np.pi)
    rpoint = rpoint_ref * width / 53
    ps = np.pi * rpoint**2
    return ps

def analyze_clusters(coordinates, labels):
    """Analyze and print details of the identified clusters."""
    coordinates['cluster'] = labels
    oligomer_sizes = coordinates['cluster'].value_counts()
    monomer_count = oligomer_sizes.get(-1, 0)  # Count of monomers
    print("Number of oligomers:", len(oligomer_sizes) - (1 if -1 in oligomer_sizes.index else 0))
    print(f"Number of monomers: {monomer_count}")
    print("\nDetails for oligomers:")
    print(f'Label : # monomers')
    for cluster_label, size in oligomer_sizes.items():
        if cluster_label != -1:  # Ignore monomers
            print(f'{cluster_label} \t {size}')

    return oligomer_sizes, monomer_count

def plot_clusters(coordinates, labels, plot_clusters, width, cluster_plot_file, color_by_cluster=0, seed=DEFAULT_SEED):
    """Visualize the clustering results."""
    newcoordinates = coordinates
    newcoordinates['cluster'] = labels
    unique_labels = set(labels)
    num_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    # Define a point size on the cluster plot to figure out protein width
    psize = calculate_point_size(coordinates, width)
    
    if plot_clusters:
        if color_by_cluster :
           hues = np.linspace(0, 1, num_clusters, endpoint=False)
           np.random.seed(seed)
           np.random.shuffle(hues)
           colors = [mcolors.hsv_to_rgb((hue, 1, 1)) for hue in hues]
           color_map = {label: colors[idx] for idx, label in enumerate(unique_labels) if label != -1}
        else :
           # Create a dictionary to hold cluster sizes
           cluster_sizes = newcoordinates['cluster'].value_counts()
           # Generate colors based on cluster sizes
           size_to_color = generate_colors_by_size(cluster_sizes)
          
        plt.figure(figsize=(10, 10))
        for label in unique_labels:
           if color_by_cluster :
               color = color_map.get(label, 'k')  # noise points/monomers are black
           else :
               if label == -1:
                  color = 'k'
               else :
                  # Get the size of the current cluster
                  size = cluster_sizes[label]
                  color = size_to_color[size]  # Get the corresponding color for the cluster size
               

           plt.scatter(newcoordinates[labels == label]['x'],
                       newcoordinates[labels == label]['y'],
                       color=color, alpha=0.7, s=psize)
        plt.title('Particle Clustering')
        plt.xlabel('X (Å)')
        plt.ylabel('Y (Å)')
        plt.grid()
        plt.axis('equal')              
        plt.grid(False)	# Remove grid lines        
        
        # Save the plot if a filename is provided
        if cluster_plot_file:
            plt.savefig(cluster_plot_file)
            print(f"Cluster plot saved as: {cluster_plot_file}")
        else:
            plt.show()

def plot_size_distribution(size_distribution, monomer_count, plot_distrib, distribution_plot_file):
    """Plot the distribution of oligomer sizes and log the results."""
    # Prepare the distribution for plotting
    size_distribution = size_distribution[size_distribution > 1]  # Only consider clusters with more than one monomer
    
    # Include monomers in the distribution
    total_monomers_distribution = size_distribution.value_counts().sort_index()
    
    # Create a new Series for monomers
    if monomer_count > 0:
        monomer_series = pd.Series({1: monomer_count})
        total_monomers_distribution = pd.concat([total_monomers_distribution, monomer_series])

    # Log the final distribution
    print("\nSize distribution of oligomers:")
    print(f'Cluster Size   # Oligomers   Total Monomers')
    for size in total_monomers_distribution.index:
        print(f'\t {size} \t {total_monomers_distribution[size]} \t {size * total_monomers_distribution[size]}')

    if plot_distrib:
        plt.figure(figsize=(10, 6))
        plt.bar(total_monomers_distribution.index, total_monomers_distribution.values, alpha=0.7)
        plt.title('Oligomer Size Distribution')
        plt.xlabel('Oligomer Size (Number of Monomers)')
        plt.ylabel('Number of Oligomers')
        plt.xticks(total_monomers_distribution.index)
        plt.grid()
        
        # Save the distribution plot if a filename is provided
        if distribution_plot_file:
            plt.savefig(distribution_plot_file, dpi=300)
            print(f"Distribution plot saved as: {distribution_plot_file}")
        else:
            plt.show()

        
def write_cluster_indices(coordinates, labels, output_file):
    """Write the indices of particles in each cluster to a text file."""
    if output_file:
        unique_labels = set(labels)
        with open(output_file, 'w') as f:
            for label in unique_labels:
                if label != -1:  # Ignore noise points
                    indices = np.where(labels == label)[0]  # Get the indices of the particles in this cluster
                    f.write(' '.join(map(str, indices)) + '\n')

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Clustering of particle within agregates using the DBSCAN algorithm ')
    parser.add_argument('--eps', type=float, default=DEFAULT_EPS, help='Maximum radius of particles to be considered in the neighborhood search.')
    parser.add_argument('--min_samples', type=int, default=DEFAULT_MIN_SAMPLES, help='Minimum size for a cluster.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input file containing particle coordinates.')
    parser.add_argument('-o', '--output_file', type=str, help='Path to the output file for cluster indices.')
    parser.add_argument('-c', '--plot_clust', action='store_true', help='Whether to plot the cluster (False without switch).')
    parser.add_argument('-d', '--plot_dist', action='store_true', help='Whether to plot the distribution of clustering (False without switch).')
    parser.add_argument('-w', '--width', type=float, default=DEFAULT_DIAM, help='Protein width for ploting cluster (in A), eg. diameter')
    parser.add_argument('--colbyclust', action='store_true', help='Each cluster is colored with distinct colors rather than being represented by size')
    parser.add_argument('-s', '--seed', type=int, default=DEFAULT_SEED, help='Seed for color plot.')
    parser.add_argument('--cpng', type=str, help='Path to save the cluster plot as an image file.')
    parser.add_argument('--dpng', type=str, help='Path to save the distribution plot as an image file.')

    # Parse the arguments
    args = parser.parse_args()
    color_by_cluster = args.colbyclust
    # Load data and perform clustering
    coordinates = load_data(args.input_file)
    distances = calculate_distances(coordinates)
    labels = perform_clustering(distances, args.eps, args.min_samples)

    # Plot clusters and analyze results
    plot_clusters(coordinates, labels, args.plot_clust, args.width, args.cpng, color_by_cluster )
    oligomer_sizes, monomer_count = analyze_clusters(coordinates, labels)
    plot_size_distribution(oligomer_sizes, monomer_count, args.plot_dist, args.dpng)
    
    # Write cluster indices to output file if specified
    write_cluster_indices(coordinates, labels, args.output_file)

if __name__ == "__main__":
    main()
