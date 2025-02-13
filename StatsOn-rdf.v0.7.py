import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Analysis of g(R) files to characterize the heterogeneity within an oligomeric aggregate.')
    parser.add_argument('-p', '--path', type=str, default='./RDF',
                        help='Path to reach the individual RDF files')
    parser.add_argument('-d', '--distrange', type=str, default='47-49 51-55 57-61 63-67 69-79 81-85',
                        help='Distance ranges to group as a single bar (blank separated fields)')
    parser.add_argument('-s', '--size', type=int, default=100,
                        help='Number of patches to include in the statistical analysis')
    parser.add_argument('-sd', '--seed', type=int, default=1970,
                        help='Seed for the random selection of patches')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    path = args.path
    size = args.size
    seed = args.seed
    
    # Convert distance ranges into a list of tuples
    dist_ranges = [tuple(map(int, range.split('-'))) for range in args.distrange.split()]
    
    #========================================================
    # Retrieve data from RDF files and tabulate them as a dataframe
    #========================================================
    	
    # Initialize an empty dataframe to store raw data and processed data
    data_brut = pd.DataFrame(columns=['dist'])
    merge_dist = pd.DataFrame(columns=['dist'])
    filtered_df = pd.DataFrame()
    
    # Initialize a list to store filtered data
    filtered_data = []
    dist_data = []

    # Retrieve the names of all .dat files in the current folder
    # Check if the specified path exists
    if not os.path.isdir(path):
       print(f"The specified path {path} does not exist.")
       exit(1)  # Exit if the path is not found

    
    # Global RDF
    fileref = os.path.join(path, "000.ref")
    # Check if the reference file exists
    if not os.path.isfile(fileref):
       print(f"Reference file {fileref} does not exist.")
       exit(1)  # Exit if the reference file is not found

    # RDFs for patches
    files = glob.glob(os.path.join(path, "*.dat"))
    # Check if any .dat files were found
    if not files:
        print(f"No .dat files found in the directory: {path}")
        exit(1)  # Exit if no .dat files are found   
        
    # The reference file will be systematically at the first data row of the dataframe
    files.insert(0, fileref)
    n = len(files)
    print(path, n)
    
    # Resize sampling to the maximum number of files if sampling > number of files
    size = min(size,n-1)

    # Random selection of indexes that will correspond to the selection of RDF files
    np.random.seed(seed)
    Alea = np.random.choice(np.arange(1,n), size=size, replace=False)
    Alea = np.insert(Alea, 0, 0)

    # Initialize an empty dataframe to record RDF data
    data_rdf = pd.DataFrame(columns=['dist'])
    # The first column of the dataframe is the reference RDF 
    try:
        df = pd.read_csv(fileref, delimiter=' ', header=0)
    except Exception as e:
        print(f"Error reading reference file {fileref}: {e}")
        exit(1)  # Exit if the reference file cannot be read

    # The other selected files for patches are then retrieved
    for n in Alea:
        file = files[n]
        try:
            df = pd.read_csv(file, delimiter=' ', header=0)
            # The first column in the .dat file contains distances
            dist = df.iloc[:, 0]
            # Check if the first column of the dataframe has already been initialized
            if len(data_rdf) == 0:
                # Initialize the first column with the distances from the .dat file
                data_rdf['dist'] = dist            
                # Check if the distances from the .dat file are identical to those in the dataframe
                if np.allclose(dist, data_rdf['dist']):
                    # If yes, add the data from the .dat file to the dataframe
                    data_rdf = pd.concat([data_rdf, df.iloc[:, 1:]], axis=1)
                else:
                    # If not, print an error message
                    print(f"Distances in {file} file do not match those of the dataframe.")
            else:
                # Check if the distances from the .dat file are identical to those in the dataframe
                if np.allclose(dist, data_rdf['dist']):
                    # If yes, add the data from the .dat file to the dataframe
                    data_rdf = pd.concat([data_rdf, df.iloc[:, 1:]], axis=1)
                else:
                    # If not, print an error message
                    print(f"Distances in {file} file do not match those of the dataframe.")
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue  # Skip to the next file if an error occurs
    
    # Add the data from the current folder to the raw data
    df_concat = data_rdf.merge(data_brut.iloc[:, 1:], left_index=True, right_index=True, how='left')
    data_brut = df_concat
    
    #================================================================================================================
    # r ratio: a tabulated conversion factor to retrieve populations from RDF (#TODO: increment assumed to be 2 Å)
    #===============================================================================================================
    for index, distance in data_rdf['dist'].items():
        try:
            # The r ratio applies the concentric sphere slice weighting to retrieve population numbers from the RDF
            r2 = float(distance) + 1 	#TODO: external slice (caution: the hard-coded integers assume a 2 Å concentric layer)
            r1 = r2 - 2			#TODO: internal slice 
            dvol = 4/3 * np.pi * ( r2 ** 3 - r1 **3) # Volume difference
            data_rdf.loc[index, 'r'] = dvol / 30.0   # Normalization to reach the number of lipids in the integral of g(r)	
        except ValueError as e:
            print(f"ValueError for distance index {index}: {e}")
            data_rdf.loc[index, 'r'] = np.nan
        except Exception as e:
            print(f"Error calculating r ratio for distance index {index}: {e}")
            data_rdf.loc[index, 'r'] = np.nan

    # Get the column names you want to update
    cols_to_update = data_rdf.columns.difference(['dist', 'r'])

    # Update the columns to convert RDF into numbers
    for col in cols_to_update:
        data_rdf[col] = data_rdf[col] * data_rdf['r']        
    del data_rdf['r']
    
    #=========================================================================
    # Sum and concatenate lines corresponding to the specified distance range
    #=========================================================================
    for dist_range in dist_ranges:
        # Filter the data based on the distance range
        filtered = data_rdf[data_rdf['dist'].between(*dist_range)]
        avdist = filtered['dist'].mean(axis=0)
        # Sum the filtered data by columns            
        filteredsum = filtered[cols_to_update].sum(axis=0)
        # Add a new index-value pair at the beginning of the Series            
        dist_series = pd.Series([avdist], index=['dist'])
        # Add the filtered data to the list of filtered data
        filtered_data.append(filteredsum)
        dist_data.append(dist_series)

    
    # Create a dataframe to store the filtered data
    dist_df = pd.DataFrame(dist_data)
    filtered_df = pd.DataFrame(filtered_data)
    
    #================================================================================
    # Extract means and standard deviations from the lines of the filtered dataframe within each distance range
    #================================================================================
    # Exclude the column at index 1 (the reference column) for statistical analysis
    filtered4stat_df = filtered_df.drop(filtered_df.columns[0], axis=1)
    cols_to_update4stat = cols_to_update.drop('RDF')

    # Calculate statistics for each row
    dist = dist_df.dist
    mean = filtered4stat_df[cols_to_update4stat].where(filtered4stat_df[cols_to_update4stat] > 0).mean(axis=1)
    std = filtered4stat_df[cols_to_update4stat].where(filtered4stat_df[cols_to_update4stat] > 0).std(axis=1)

    # Create a dataframe to store the statistics
    stat_df = pd.DataFrame({'Dist': dist, 'mean': mean, 'std': std})              

    # Print the final statistics
    print(stat_df)
    
    ##===================================
    # Graph production using Matplotlib
    #====================================
    # Graph configuration
    fig, axs = plt.subplots(nrows=1, figsize=(10, 4.9))
    ax = plt.gca()  # Get current axes
    # Remove the frames
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.yaxis.tick_right()  # Move y-ticks to the right

    # Adjust the position of labels on the x-axis
    axs.set_xticks([j*5 for j in stat_df['Dist']])
    axs.set_ylim(0, 3.0)

    # The bars are processed one by one in a loop over the distances
    for j in range(len(stat_df['Dist'])):
        axs.bar(j*3, filtered_df['RDF'][j], width=1.0, color='#b81b68') # Red 'color blind' for global RDFs
        axs.bar(j*3-0.5, stat_df['mean'][j], width=1.0, color='#7898e7') # Blue 'color blind' for patch RDFs
        # Variability of the patches
        axs.errorbar(j*3-0.5, stat_df['mean'][j], yerr=stat_df['std'][j], fmt='none', capsize=3, capthick=1, elinewidth=1, color='black')
        
        # Overlay individual values as dots
        data = filtered_df[cols_to_update4stat].where(filtered_df[cols_to_update4stat] > 0).iloc[j]
        xdata = np.full_like(data, j*3-0.5)
        axs.scatter(xdata, data, color='gray', label='Individual Values')

    # Display the graph
    plt.show()
