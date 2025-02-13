#!/usr/bin/python
# -*- coding:LATIN-1 -*-
# Importer les modules nécessaires
import argparse
#import cupy as cp
import numpy as np
import pandas as pd
import MDAnalysis as mda
import warnings
import sys
import itertools
import multiprocessing
import re

if sys.hexversion < 0x03000000:
    mapper= itertools.imap # 2.4 ² Python < 3
else:
    mapper= map # Python ³ 3


def loadgro(files, file_numbers):
    tous_les_atomes = []
    syst = mda.Universe(files[0])
    box = syst.dimensions
    boxes = []
    dcyl = dcyzone 
    
    # Positionner les protéine suivant les coordonnées x y des données
    for i in range(len(files)):
        pnum = file_numbers[i]
        indx = pnum - 1
        print(f"Loading and positioning Prot#{pnum}")
        #Chargement du nieme fichier gro
        syst= mda.Universe(files[i])
        #Attribution d'une identification de segment
        Attrib= ['p{:03d}'.format(pnum)]
        syst.add_TopologyAttr('segid', Attrib)
        # Obtenir la dimension de la boÃ®te des conditions périodiques
        boxnew = syst.dimensions
        #Restreindre la sélection à un cylindre autour du COG de la protéine
        max_attempts = 1
        attempt = 0
        while attempt <= max_attempts :
                try :
                        sel = "same resid as cyzone "+ str(dcyl) + f" {boxnew[2]/2 } -{boxnew[2]/2} protein"
                        syst=syst.select_atoms(sel, periodic = 'False')
                        break
                except NotImplementedError as e :
                        print("The default cyzone radius : " + str(dcyzone) + " is replaced by the maximum value corresponding to the x PBC ", boxnew[1])
                        dcyl = boxnew[1]
                        attempt = attempt + 1
                except Exception as e :
                        print("Unrecoverable error, exiting")
                        sys.exit()
        
        ##Positionnement de la protéine suivant les coordonnées AFM
        #Récupération des info du fichier de positions
        x, y, z = positions[indx][0:3]
        # Appliquer la translation
        syst.atoms.translate([x, y, z])
        # Nouvelle boite en prenant le maximum des dimensions dans chaque direction et la translation
        box = np.maximum(box, boxnew + [x, y, z,0,0,0])  
        # Export des coordonnees vers l'univers global
        tous_les_atomes.append(syst.atoms)
    boxes.extend(box)
    return tous_les_atomes, boxes

#Calcule l'angle entre l'axe x du systeme et une face protéique de reférence
def calcangle(univ, numatom):
    prot = univ
    # Calculer le vecteur entre le centre géométrique et l'atome de référence
    v = prot[num_atom].position[0:2]
    #The Arctan2 function directly give the angle in respect with the  x axis
    angle = np.arctan2(v[1], v[0])
    if angle < 0:
        angle += 2 * np.pi
    return angle

def rotsys(univ, num_atom, theta):
    #Could apply either on a full syst or protein only selection, depending of the content of univ object
    #This is to assure that a subselection is only made of the protein content of the system
    prot = univ.select_atoms("protein", periodic = 'False')
    cg = prot.center_of_geometry()
    #Locate the prot/sys in the origine to prepare the rotation along the protein z axis
    univ.atoms.translate(-cg)    
    #Calcule l'angle entre l'axe x du systeme et une face protéique de reférence
    angle=calcangle(prot, num_atom)
    # Convertir theta de rotation en radians
    theta = np.radians(theta)  
    myrot = theta - angle
    
    # Créer une matrice de rotation autour de l'axe z
    R = np.array([[np.cos(myrot), -np.sin(myrot), 0],
                  [np.sin(myrot), np.cos(myrot), 0],
                  [0, 0, 1]])
    # Apply rotation to prot/sys
    univ.rotate(R)
    univ.atoms.translate(cg)  
    return univ

# Create a dataframe 
def mkdf(name):
    name = pd.DataFrame({
    "protein_number": [],
    "angle_value": [],
    "clash_number": [],
    "Segids": [],
     "Segids2": []
})
    return name

#Filter a bestrot dataframe to select best rotation exhibiting minimal claches
def filtrot(bestrot):
    grouped_df = bestrot.groupby("protein_number")
    min_idx = grouped_df["clash_number"].idxmin()
    filtered_rot = bestrot.loc[min_idx]   
    return  filtered_rot   

#Unclashing by translation #Dev : not working !
def unclashp(Attrib, univ, lsegids):
    #Could apply either on a full syst or protein only selection, depending of the content of univ object
    #This is to assure that a subselection is only made of the protein content of the system
    prot = univ.select_atoms("protein", periodic = 'False')
    prot1 =  prot.select_atoms("segid "  + "".join(Attrib), periodic = 'False')
    cg1 = prot1.center_of_geometry()
    nclash = 1
    for i, seg in enumerate(lsegids):
        sel = "segid " + "".join(seg) + " and around " + str(dist)  +" segid " + "".join(Attrib)
        while nclash > 0  :
            clash = prot.select_atoms(sel, periodic = 'False')            
            cg2 = clash.center_of_geometry()
            trans=cg1-cg2
            prot1.translate(trans)
            print(len(clash), trans)            

     
def findneigbor_mult(shunksegid):
    distnb = dcyzone * 2 #twice dcyzone in loadgro (a maximum corresponding to PBC may limit the choice)
    mysegids = shunksegid
    minunivers = univers.select_atoms('protein', periodic = 'False') #A mettre avant l'appel de la fonction ?
    index=0
    bestrot = mkdf("bestrot")
    
    #Loop on proteins

    for Attrib in sorted(mysegids) :
        print("Searching clashing neighbor for ", Attrib)
        prot_num = int(Attrib.replace("p", ""))
        indx = prot_num -1
        #Select atoms for the whole segname and for its protein only
        mysegid = " segid "+ "".join(Attrib)
        #Selection for chashes of all other segids over the current one
        #For CG
        protatom = " name BB and "
        #For All at
        #protatom = " name CA CB and "
        selclashp = protatom + " not" + "".join(mysegid) + " and around " + str(dist) + protatom + "".join(mysegid)
        selclashp2 = protatom + " not " + "".join(mysegid) + " and around " + str(distnb) + protatom + "".join(mysegid)

        clashp = minunivers.select_atoms(selclashp, periodic = 'False')
        clashp2 = minunivers.select_atoms(selclashp2, periodic = 'False')
        nbclash = len(clashp)
        segidlist = set(clashp.segids)
        segidlist2 = set(clashp2.segids)
        segidlist.add(Attrib)
        segidlist2.add(Attrib)
        # The results of the attempt on the dataframe
        bestrot.loc[index, ["protein_number","angle_value","clash_number", "Segids", "Segids2"]] = [prot_num, theta_list[indx], nbclash, segidlist, segidlist2 ]
        index = index + 1
        bestrot = filtrot(bestrot)
    return bestrot
 
def minclashmulti(arg):
    sel = "segid " + " ".join(arg)
    mysegids = list(arg)
    prot_num = [int(re.search(r"p(\d+)", s).group(1)) for s in mysegids]
    indx = [x - 1 for x in prot_num]
    minunivers = univers.select_atoms(sel, periodic = 'False')
    index = 0

    # #awidth = 360/nbrot
    awidth = 1    
    #Loop on proteins
    for Attrib in mysegids :
        searchrot = mkdf('searchrot')
        print("Minimizing Clashes for ", Attrib)
        #To retriev the angle coresponding to the index of prot indexed as j
        j = 0
        k = indx[j]
        angle=theta_list[k]
        theta=angle
        #Select atoms for the whole segname and for its protein only
        mysegid = " segid "+ "".join(Attrib)
        segsyst = minunivers.select_atoms(mysegid, periodic = 'False')
        psyst = segsyst.select_atoms("protein", periodic = 'False')
        #Selection for chashes of all other segids over the current one
        selclashp = "name BB and not" + "".join(mysegid) + " and around " + str(dist) + " name BB and" + "".join(mysegid)
        #Rotate by awidth and look for clashes
        
        rotmax = int(360/awidth)
        nrot = 0
        while nrot < rotmax :
            #Rotation about z
            psyst=rotsys(psyst, num_atom, theta)
            clashp = minunivers.select_atoms(selclashp, periodic = 'False')            
            nbclash = len(clashp)
            #Building a set of clashing prot
            segidlist = set(clashp.segids)
            segidlist.add(Attrib)
            #Replace the proteine in it's original location
            psyst=rotsys(psyst, num_atom, -theta)
            # The results of the attempt on the dataframe
            searchrot.loc[index, ["protein_number","angle_value","clash_number", "Segids"]] = [Attrib, theta%360, nbclash, segidlist ]
            index = index + 1
            if nbclash == 0 : break
            #After the first rotation the next will be additive by awidth
            theta = theta + awidth
            nrot+=1   
        #Puts the segid sys to the best found rotation
        segsyst=rotsys(segsyst, num_atom, theta)
    searchrot = filtrot(searchrot)
    return searchrot
   

def chunking(lst, n):
    #Yield successive n-sized chunks from lst#
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
def split_list(lst, n):
    """Splits a list into sublists of size n."""
    return [lst[i:i+n] for i in range(0, len(lst), n)]

#Provide a list of nonredundant cluster list based on the proximity assesment adressed by findneigbor_mult()
def hierarchical_fusion(clusters):
    # Initialize the final clusters list
    final_clusters = []

    # Sort the clusters by length in descending order
    clusters.sort(key=len, reverse=True)

    for cluster in clusters:
        # Convert each cluster to a set for efficient comparison
        cluster_set = set(cluster)

        # Check if the current cluster overlaps with any final cluster
        overlap = False
        for final_cluster in final_clusters:
            if cluster_set & final_cluster:  # If there's an intersection
                final_cluster |= cluster_set  # Merge the two clusters
                overlap = True
                break

        # If it doesn't overlap with any final cluster, add it to the final clusters
        if not overlap:
            final_clusters.append(cluster_set)

    # Convert the sets back to lists
    final_clusters = [list(cluster) for cluster in final_clusters]
    

    return final_clusters

def balance_clusters2(clusters, num_clusters):
    # Calculate the total number of elements and the ideal size of each cluster
    total_elements = sum(len(cluster) for cluster in clusters)
    ideal_size = total_elements // num_clusters
    max_size = int(1.5 * ideal_size)
    min_size =  int(ideal_size // 1.5)
    min_clust = num_clusters -1
    max_clust = num_clusters +1
    print("total_elements : ",total_elements," ideal elt by clust : ",ideal_size, "max_elt by clust : ",max_size,"min_elt by clust : ",min_size)
    print("min cluster number : ",min_clust, "max clust number : ",max_clust)

    # Convert each cluster to a list and sort them by size in descending order
    clusters = sorted([list(cluster) for cluster in clusters], key=len, reverse=True)


    # Check if any cluster is larger than the maximum size or smaller than the minimum size
    for i in range(len(clusters)-1):	#To allow fusion with the true last cluster
    	#Trop grand : division en deux du plus grand
        while len(clusters[i]) > max_size:
            print("Cluster size : ", len(clusters[i]))
            # If a cluster is too large, split it into two
            half = ideal_size
            clusters.append(clusters[i][:half])
            clusters[i] = clusters[i][half:]
        #Trop petit : fusion avec un autre  plus petit     
        while len(clusters[i]) < min_size and len(clusters) > 1:
            print("Cluster size : ", len(clusters[i]))
            if len(clusters[i]) < 2 : break
            # 
            other_clusters = clusters[:i] + clusters[i+1:]
            smallest_other_cluster = min(other_clusters, key=len)
            clusters.remove(smallest_other_cluster)
            clusters[i] += smallest_other_cluster

    # Test the last cluster for a posible fusion with the previous one
    clusters.sort(key=len, reverse=True)
    j = len(clusters) -1
    if len(clusters[j]) < min_size :
	    other_clusters = clusters[:j] + clusters[j-1:]
	    smallest_other_cluster = min(other_clusters, key=len)
	    clusters.remove(smallest_other_cluster)
	    clusters[i-1] += smallest_other_cluster
 

    # Sort the clusters by size in descending order
    clusters.sort(key=len, reverse=True)

    # While the number of clusters is not within the desired range
    while not (min_clust <= len(clusters) <= max_clust):
        if len(clusters) < 3 : break
        print("Num cluster : ", len(clusters), min_clust, max_clust )
        # If there are too many clusters
        if len(clusters) > max_clust :
            i=0
            #Zap the cluster if aready in ideal size
            while len(clusters[i])  == ideal_size :
                i = i+1
            # Pop the largest cluster not ideally sized  
            largest_cluster = clusters.pop(i)                
            # If the largest cluster is bigger than the ideal size, split it into two
            if len(largest_cluster) > ideal_size:
                #half = len(largest_cluster) // 2
                clusters.append(largest_cluster[:ideal_size])
                remain_cluster = largest_cluster[ideal_size:]
            else:
                # Otherwise, find the smallest cluster and merge it with the largest
                remain_cluster = largest_cluster
            #Assume cluster fusion to be sure that the len(cluster) lower    
            smallest_cluster = min(clusters, key=len)
            clusters.remove(smallest_cluster)
            clusters.append(remain_cluster + smallest_cluster)
        else:  # If there are too few clusters
            if len(clusters) < 3 : break
            # Pop the  largest cluster if not ideally sized
            i=0
            while len(clusters[i])  == ideal_size :
                i = i+1                
            largest_cluster = clusters.pop(i)  
            # And break it in two parts              
            clusters.append(largest_cluster[:ideal_size])
            clusters.append(largest_cluster[ideal_size:])
                         
    return clusters


#Constraints the cluster list on the number of cluster to define (based on the cpu to use)
def balance_clusters(clusters, num_clusters):
    # Calculate the total number of elements and the ideal size of each cluster
    total_elements = sum(len(cluster) for cluster in clusters)
    ideal_size = total_elements // num_clusters
    max_size = int(1.5 * ideal_size)
    min_size =  int(ideal_size // 1.5)
    min_clust = num_clusters -1
    max_clust = num_clusters +1
    print(" Total elements : ",total_elements,"\n Ideal elt number by clust : ",ideal_size, "\n Max targeted elt by clust : ",max_size,"\n Min targeted elt by clust : ",min_size)
    print("Min cluster number : ",min_clust, "\n Max clust number : ",max_clust)

    # Convert each cluster to a list and sort them by size in descending order
    clusters = sorted([list(cluster) for cluster in clusters], key=len, reverse=True)

    # While the number of clusters is not within the desired range
    while not (min_clust <= len(clusters) <= max_clust):
        if len(clusters) < 3 : break
        print("Num cluster : ", len(clusters), min_clust, max_clust )
        # If there are too many clusters
        if len(clusters) > max_clust :
            # Pop the it largest cluster if not ideally sized
            i=0
            if len(clusters[i])  == ideal_size :
                i = i+1
                
            largest_cluster = clusters.pop(i)                
            # If the largest cluster is bigger than the ideal size, split it into two
            if len(largest_cluster) > ideal_size:
                #half = len(largest_cluster) // 2
                clusters.append(largest_cluster[:ideal_size])
                remain_cluster = largest_cluster[ideal_size:]
            else:
                # Otherwise, find the smallest cluster and merge it with the largest
                remain_cluster = largest_cluster
            #Assume cluster fusion to be sure that the len(cluster) lower    
            smallest_cluster = min(clusters, key=len)
            clusters.remove(smallest_cluster)
            clusters.append(remain_cluster + smallest_cluster)
        else:  # If there are too few clusters
            # Find the two smallest clusters and merge them
            if len(clusters) < 3 : break
            # Pop the it largest cluster if not ideally sized
            i=0
            if len(clusters[i])  == ideal_size :
                i = i+1                
            largest_cluster = clusters.pop(i)                
            clusters.append(largest_cluster[:ideal_size])
            clusters.append(largest_cluster[ideal_size:])
       
        # Sort the clusters by size in descending order
        clusters.sort(key=len, reverse=True)

    # Check if any cluster is larger than the maximum size or smaller than the minimum size
    for i in range(len(clusters)-1):	#To allow fusion with the true last cluster
    	#Trop grand : division en deux du plus grand
        while len(clusters[i]) > max_size:
            print("Cluster size : ", len(clusters[i]))
            # If a cluster is too large, split it into two
            half = len(clusters[i]) // 2
            clusters.append(clusters[i][:half])
            clusters[i] = clusters[i][half:]
        #Trop petit : fusion avec un autre  plus petit     
        while len(clusters[i]) < min_size and len(clusters) > 1:
            print("Cluster size : ", len(clusters[i]))
            if len(clusters[i]) < 2 : break
            # 
            other_clusters = clusters[:i] + clusters[i+1:]
            smallest_other_cluster = min(other_clusters, key=len)
            clusters.remove(smallest_other_cluster)
            clusters[i] += smallest_other_cluster

    # Test the last cluster for a posible fusion with the previous one
    clusters.sort(key=len, reverse=True)
    j = len(clusters) -1
    if len(clusters[j]) < min_size :
	    other_clusters = clusters[:j] + clusters[j-1:]
	    smallest_other_cluster = min(other_clusters, key=len)
	    clusters.remove(smallest_other_cluster)
	    clusters[i-1] += smallest_other_cluster
                   
    return clusters
    
def rm_dupli(nested_list):
    seen = set()
    new_list = []
    for sublist in nested_list:
        new_sublist = []
        for item in sublist:
            if item not in seen:
                seen.add(item)
                new_sublist.append(item)
        new_list.append(new_sublist)
    return new_list

def renumber_nonprot():
    #Renumérotation des résidus
    forprotatsize = univers.select_atoms("protein and segid p001", periodic = 'False')
    protatsize = len(forprotatsize.atoms)
    last_presid=nb_proteines * protatsize
    print("Residue renumbering")
    residues = univers.select_atoms('not protein', periodic = 'False')
    for i, residue in enumerate(residues.residues):
      residue.resid = last_presid + i + 1        
    last_resid = last_presid + i + 1
    
    return last_resid 

def delclash(args):
    univers, clash_res = args
    sel=("not protein and resid " + " ".join(str(res) for res in clash_res))  
    clash=univers.select_atoms(sel, periodic='False')   
    print("Deleting ", len(clash), " clashing atoms")
    univers.atoms -= clash
    mylen=len(univers.atoms)
    
    return  mylen

def findclashes(args):
    univers_atoms, filenumbers = args
    pindex = [x - 1 for x in filenumbers]
    clashes= univers_atoms.select_atoms("resid 0", periodic = 'False')
    for i in pindex :
        #Les clash sont recherchés entre le segid nouvellement chargée et tous les segids précédent (immunité de la prot + lipides shell 1). 
        mynewsegid = "segid " + "".join(lseg[i])
        print("Finding clashing proteins with segid ", lseg[i])
        itsanulus = " and not (protein or same resid as around 5 protein) "
        resclashothersegid =  "and same resid as (around " + str(dist) + " not segid " + " ".join(lseg[i:]) + ")"
        searchLO = "not protein and " + mynewsegid + itsanulus + resclashothersegid
        
        # clash de la lipide shell1 nouvelle avec d'autre lipides -> c'est l'autre lipides qui doit être enlevé
        otherlip=" same resid as (not segid " + "".join(lseg[i])
        myliparround=" and around " + str(dist-1.5) + " (segid " + "".join(lseg[i]) + " and not protein"
        arroundmyprot=" and around 5 (protein and segid " + "".join(lseg[i]) + ")))"
        searchLL = "not protein and " + otherlip + myliparround + arroundmyprot   
         # Liste des recherches
        searches = [searchLO, searchLL]           
        # Liste des clashes
        # Boucle pour identifier les atomes en clash
        for search in searches:
            clash = univers_atoms.select_atoms(search, periodic='False')
            clashes = clashes | clash       
        # Boucle pour lister les résidus en clashes
        m2del = clashes.residues.resids
    return m2del

def clust2unclash(arg):
    segid1, segid2 = arg
    sel = "segid " + " ".join(segid1) + " and around " + "".join(str(dist))+ " segid " + " ".join(segid2)
    clash = univers.select_atoms(sel, periodic = 'False')
    print("Found ", len(clash), " clashing non-protein atoms")
    m2del = clash.residues.resids
    mergesegid = segid1 + segid2
    
    return  mergesegid, m2del

def lipinpore2del(shunksegid):
    mysegid = shunksegid
    mysearch2=f"same resid as (resname CHOL POPC POPE and cyzone 12.5 {box[2]/2} -{box[2]/2} (protein and segid " + "".join(mysegid) + "))"
    clash=univers.select_atoms(mysearch2, periodic='False')
    m2del = clash.residues.resids
    print(len(m2del), " lipid residues next to the pore of "+ "".join(mysegid) + " to delete")

    return m2del

if __name__ == '__main__':

    # suppress some MDAnalysis warnings when writing PDB files
    warnings.filterwarnings('ignore')
    
    # Créer un objet ArgumentParser
    parser = argparse.ArgumentParser(description='Maper des fichiers de coordonnée Gromacs sur des coordonnées d\'AFM')
    
    # Ajouter des arguments
    parser.add_argument('-i', '--prefix', type=str, default='sys', help='Préfixe du nom des fichiers de coordonnées Ã  concaténer (défaut: "sys")')
    parser.add_argument('-p', '--path', type=str, default='/home/jpierre/Documents/Manips/VDAC/AFM2CG/Assemblies/Ms2Grofiles', help='Path to the input "gro" files (défaut: see code)')
    parser.add_argument('-n', '--n_prot', type=int, default=253, help='Nombre de protéines a concatener (défaut: 5)')
    parser.add_argument('-t', '--rand', action="store_true" , help='True = Generate a random orientation else use the input file')
    parser.add_argument('-nm', '--nomin', action="store_false" , help='True = Search an optimumum rientation to minimize inter protein clashes')    
    parser.add_argument('-r', '--ref', type=int, default=158, help='index de l\'atome servant a spécifier l\'orientation de référence de la protéine (défaut: 159 - E72)')
    parser.add_argument('-d', '--dist', type=float, default=3.3, help='Distance pour la détermination des clash steriques (default = 3.3)')
    parser.add_argument('-c', '--dcyl', type=float, default=50.0, help='Distance for cylinder radius to concerve around each protéin (default = 50 A)')
    parser.add_argument('-a', '--afmcoor', type=str, default='coord.dat',  help='nom du fichier de coordonnée et d\'orientation des proteines mappée sur le cliché AFM (défaut: coord.dat)')
    parser.add_argument('-o', '--output', type=str, default='catCG2AFM.gro', help='Fichiers de coordonnées de sortie ("gro" ou "pdb", défaut: "catCG2AFM.gro")')
    
    # Parse les arguments
    args = parser.parse_args()
    
    # Définir les paramètres d'entrée
    nb_proteines = args.n_prot # Nombre de protéines Ã  positionner (first index at 0)
    rot = args.rand
    minprot=args.nomin
    #minprot=True
    
    grofilespath = args.path + "/"
    prefix = args.prefix
    # Liste des fichiers de coordonn�es au format gromacs
    fichiers_gro = [grofilespath + prefix + '{:04d}.gro'.format(i) for i in range(1, nb_proteines+1)]
    print(fichiers_gro)
    
    num_atom = args.ref # Numéro d'atome qui définit l'axe de référence avec le centre géométrique de la protéine
    dist = args.dist
    dcyzone = args.dcyl
    
    fichier_dat = args.afmcoor # Fichier de donnes qui contient la position et l'angle de rotation de chaque protéine
    
    # Créer un univers vide qui contiendra les protéines positionnées
    univers = mda.Universe.empty(n_atoms=0, n_residues=0, atom_resindex=None, trajectory=True)
    
    # Lire le fichier de donnees et stocker les orientations et les positions et dans des listes et un tableau numpy
    positions = np.loadtxt(fichier_dat)
    prot_list = list(range(1,nb_proteines+1))
    theta_list=positions[:, 3].tolist()
    
    #If the random switch is set replace the list of predifined angle by a random one
    if rot :
        theta_list = np.random.randint(1, 360, len(theta_list))
        theta_list = theta_list.tolist()
    
    #Load the gro coordinates using a multicpu parser
    num_cpus = multiprocessing.cpu_count() - 2
    # Create a pool of worker processor (concider processes=None for automatic num_cpu assignement)
    pool = multiprocessing.Pool(processes=None)
    # Combine the files list and the file number list into a single list of tuples
    mychunks = list(zip(chunking(fichiers_gro, num_cpus), chunking(prot_list, num_cpus)))
    # Positionning proteins on the AFM coordinate using parallel processing
    results = pool.starmap(loadgro, mychunks)
    # Close the pool and wait for all worker processes to finish
    pool.close()
    pool.join()         
    # Get the results of parrallel processing (tuples, containing indivisual processor work)    
    listatoms, boxes = zip(*results)
    
    #Reconstitute universe from atoms
    print("Merging ", nb_proteines, "proteins in a single coordinate system")
    univers = mda.Merge(*[element for sublist in listatoms for element in sublist])
    lseg=univers.segments.segids

    print("Merge OK")
    box = [0,0,0,0,0,0]
    for boxnew in boxes:
        box = np.maximum(box, boxnew)
    univers.dimensions = box
 
    #Backup residue number of proteins
    resprotlist = []
    residprot = univers.select_atoms('protein', periodic = 'False')
    for i, residue in enumerate(residprot.residues):
      resprotlist.append(residue)

    #Renumbering non protein atoms
    last_resid = renumber_nonprot()
        
    #Seecking favorable relative orientation  
    if minprot :
        num_cpus = multiprocessing.cpu_count() - 5
        print("Looking for optimal protein orientations")
        prot_atoms = []
        prot=univers.select_atoms("protein", periodic = 'False')
        prot_atoms.append(prot.atoms)
        protuniv = mda.Merge(prot)
        segids = set(protuniv.atoms.segids)
        segids = sorted(list(segids))
        print("Finding a list of neighbors for each Protein")
        segidchunk = chunking(segids,num_cpus)
        
        #Retrive a list of conflicting prots
        # Create a pool of worker processes
        pool = multiprocessing.Pool(processes=None)
        bestrot = pool.map(findneigbor_mult, segidchunk)
        pool.close()
        pool.join()
        combined_df = pd.concat(bestrot)
        bestrot = combined_df.reset_index(drop=True)
        
        #Extract list of conflicting prots
        seg_list = bestrot["Segids"].tolist()
        clashlist = [list(seg_list[i]) for i in range(len(seg_list)) if len(seg_list[i]) > 1]
        clashlist = sorted({tuple(sorted(item)) for item in clashlist})
        
        #Multi processing of clash resolution
        if len(clashlist) >= multiprocessing.cpu_count() - 2 :
            num_cpus = multiprocessing.cpu_count() - 2  
        else :
            num_cpus = len(clashlist) + 1            
        pool = multiprocessing.Pool(processes=None)      
        nbestrot=pool.map(minclashmulti, clashlist, chunksize=1)
        pool.close()
        pool.join()        
        
        #Save the bestrot table
        bestrot.to_csv('bestrot2.csv', index=True)
        
        #Upgrade the orientation list with the optimised orientations
        print("Applying optimal orientations")             
        prot_list = bestrot["protein_number"].tolist()
        prot_list = pd.Series(prot_list)
        prot_list = prot_list.astype(int)
        theta_list = bestrot["angle_value"].tolist()
        
    ##Impose orientation (read from input or optimised)
    for i in range(nb_proteines):
        ## NUméro de la première protéine
        pnum = prot_list[i]
        ptheta = theta_list[i]
        print(f"Prot#{pnum} final orientation :" , ptheta) 
        Attrib = ['p{:03d}'.format(pnum)]
        selrot = "segid " + "".join(Attrib)
        sysrot = univers.select_atoms(selrot, periodic = 'False')       
        sysrot=rotsys(sysrot, num_atom, ptheta)

    print('Retrieving clashes between neighbor segids')
    #Extract list of conflicting prots
    seg_list2 = bestrot["Segids2"].tolist()
    clashlist2 = [list(seg_list2[i]) for i in range(len(seg_list2)) if len(seg_list2[i]) > 1]
    clashlist2 = sorted({tuple(sorted(item)) for item in clashlist2})
       
    print('Initiating hierachic clustering for finding clashes')
    sel=[]
    indx=[]
    mysegids=[]
    prot_num=[]
    args4fc=[]
    alluniv_atoms=[]
    i=0

    # Sorting of proteins over clashing neighbor and supression of duplicate entries
    #ini_clusters = hierarchical_fusion(clashlist2)
    #ini2_clusters = rm_dupli(ini_clusters)
    ini2_clusters = rm_dupli(clashlist2)
    #Preparing multiprocessing
    maxcpu = multiprocessing.cpu_count()
    memblowguard = 0.25 # To avoid crash over memory lack
    Availcpu = int(maxcpu // (1 + memblowguard))
    
    tot_elt = len([elt for segs in ini2_clusters for elt in segs])
    wkload = tot_elt // Availcpu
    eltbyproc = 10
    if wkload >= 10 :
        num_clusters = 20
        num_cpus = Availcpu
    elif 5 <= wkload < 10 :
        num_clusters = 10
        num_cpus = Availcpu
    elif 1 <= wkload < 5 :
        num_clusters = 10
        num_cpus = num_clusters
    else:
        num_clusters = 1
        num_cpus = 1 
        
    # re-distribution of work over num cluster 
    #final_clusters = balance_clusters(ini2_clusters, num_clusters)   #Test the optimum number of CPU for the size of the cluster
    final_clusters = ini2_clusters
    #Arguments build up for multiprocessing processing of clash resolution
    #using proximity based clusters
    for seg in final_clusters :
        sel.append("segid " + " ".join(seg))
        mysegids.append(list(seg))
        prot_num.append([int(re.search(r"p(\d+)", s).group(1)) for s in mysegids[i]])
        indx.append([x - 1 for x in prot_num[i]])
        minunivers_atoms = univers.select_atoms(sel[i], periodic = 'False')
        alluniv_atoms.append(minunivers_atoms)
        args4fc.append((alluniv_atoms[i],prot_num[i]))
        i = i +1
        
    # Launch the processing of each chunk of protein numbers in parallel
    pool = multiprocessing.Pool(num_cpus)          
    m2del = pool.map(findclashes, args4fc)
    # Close the pool and wait for all worker processes to finish
    pool.close()
    pool.join()

    # Get the results of parrallel findclash
    list2del = set([])
    if len(m2del) > 0:
        reslist = np.concatenate(m2del).tolist()
        reslist = set(reslist)
        list2del = list2del.union(reslist)


    #Merge protein selection based on clusters
    clust =  final_clusters
    while len(clust) > 1 :
        newclust=[]
        chunks=[]
        #Cluster are merged by pair till it remain a single cluster
        for i  in range(0, len(clust) - 1, 2) :
            chunks.append((clust[i],clust[i+1]))  
        pool = multiprocessing.Pool(processes=None)
        results = pool.map(clust2unclash,chunks)
        pool.close()
        pool.join()
        newclust, m2del = zip(*results)
        print(len(newclust))
        clust = newclust
        if len(m2del) > 0 :
            reslist = np.concatenate(m2del).tolist()
            reslist = set(reslist)
            list2del = list2del.union(reslist)
    

    #Searching lipid, water and ions clashing with proteins   (on a single proc)    
    sel=("not protein and same resid as around " + "".join(str(dist))+ " protein")
    clash=univers.select_atoms(sel, periodic='False')
    m2del = clash.residues.resids
    reslist = set(m2del)
    list2del = list2del.union(reslist)
   
    
    # Searching Lipids headen in protein pores 
    pool = multiprocessing.Pool(processes=None)
    m2del = pool.map(lipinpore2del,seg_list)
    pool.close()
    pool.join()  
    #m2del=rm_dupli(m2del)
    if len(m2del) > 0 :
        reslist = np.concatenate(m2del).tolist()
        reslist = set(reslist)
        list2del = list2del.union(reslist)

    list2del = sorted(list2del)
    print("Deleting ", len(list2del), " molecules")
    
    delclash((univers,list2del))
    
    #Renumérotation des résidus non proteiques
    last_resid = renumber_nonprot()
    
    #Restauration des resids prot (Really require ?)
    residprot = univers.select_atoms('protein', periodic = 'False')
    for i, residue in enumerate(residprot.residues):
      residue.resid = resprotlist[i].resid
   
    #Write final outputs
    print("Writing gro output :", args.output)
    univers.atoms.write(args.output)
    
    mypdb = "Prot" + str(nb_proteines) + ".pdb"
    print("Writing pdb output :", mypdb)
    univers.atoms.write(mypdb)  
