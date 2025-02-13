# Article MolBiol 2025
## Lipids modulate the organization of VDAC1, the gatekeeper of mitochondria Mol. Biol 2025

![CG Image](./IMG/myimg.jpg)

V1 - preprint in [BioRXive](https://www.biorxiv.org/content/10.1101/2024.06.26.597124v1)

VDACs, the most abundant proteins in the outer mitochondrial membrane (MOM), are crucial for mitochondrial physiology. VDAC regulate metabolite and ion exchange, modulate calcium homeostasis, and play roles in numerous cellular events such as apoptosis, mitochondrial DNA (mtDNA) release, and different diseases. Mitochondrial function is closely tied to VDAC oligomerization, influencing key processes like mtDNA release and apoptosis, but the molecular drivers of this oligomerization remain unclear. In this study, we investigate the effects of three major MOM lipids on VDAC assemblies using atomic force microscopy and molecular dynamics simulations. Our results show that phosphatidylethanolamine and cholesterol regulate VDAC assembly, with the formation of stable lipid-protein clusters of various size and compaction. Deviations from physiological lipid content disrupted native-like VDAC assemblies, highlighting the importance of lipid environment in VDAC organization. These findings underscore how lipid heterogeneity and changes in membranes influence VDAC function.

**AFM2CG**  : a python script to graft CG (MArtini) molecular dynamics simulations coordinates on picked Atomic Forces Microscopy coordinates

V 0.97p (beta)
- The AFM coordinates are provide as a simple coord.dat, csv (blank or tab separated fields) containing x y z rot coordinate comming from AFM image peacking. Coordinates in &angst; rot is a relative orientation (deg) relative to a reference surface (see *infra*)
*Sample*
```   192.2	343.7	0.0	100
    205.0	403.3	0.0	229
    225.5	476.0	0.0	258
    238.3	377.4	0.0	21
    284.4	512.3	0.0	257
    289.5	439.7	0.0	34
    294.6	351.5	0.0	200
```
- A dedicated directory, referencend to the "path" argument contains Coarse Grained MARTINI cordinates in a GRO frormat (Gromacs).
- Other parameter defines in the argument parser :
 
```
   #Arguments
    parser.add_argument('-i', '--prefix', type=str, default='sys', help='Préfixe du nom des fichiers de coordonnées Ã  concaténer (défaut: "sys")')
    parser.add_argument('-p', '--path', type=str, default='./CGGrofiles', help='Path to the input "gro" files (défaut: see code)')
    parser.add_argument('-n', '--n_prot', type=int, default=5, help='Nombre de protéines a concatener (défaut: 5)')
    parser.add_argument('-t', '--rand', action="store_true" , help='True = Generate a random orientation else use data in the input file')
    parser.add_argument('-nm', '--nomin', action="store_false" , help='True = Search an optimumum rientation to minimize inter protein clashes')    
    parser.add_argument('-r', '--ref', type=int, default=158, help='Atome reference (index) on the CG coordinate to define the reference face for orientation (0 deg) (défaut: 158 - E72)')
    parser.add_argument('-d', '--dist', type=float, default=3.3, help='Distance for clashes consideration (default = 3.3)')
    parser.add_argument('-c', '--dcyl', type=float, default=50.0, help='Distance for cylinder radius to concerve lipids around each protéin (default = 50 A)')
    parser.add_argument('-a', '--afmcoor', type=str, default='coord.dat',  help='nom du fichier de coordonnée et d\'orientation des proteines mappée sur le cliché AFM (défaut: coord.dat)')
    parser.add_argument('-o', '--output', type=str, default='catCG2AFM.gro', help='Fichiers de coordonnées de sortie ("gro" ou "pdb", défaut: "catCG2AFM.gro")')
```
*Features*
- Use python multi-threads provided by the "multiprocessing" library
- Try a hierachical clutering analysis to aptimize processor loads (bigs still present)
- Suppress lipid overlaps (either to proteins or to previously positionned lipids)
- Perform a very basic protein-protein optimization procedure to limite clashes by protein rotation

*Bugs an todo things*
- The hierachical clusterization procedure sometime behave poorly giving rise to unresolved situation vith no convergence
- Allows small translations in the peeking error range to try to resolve additionnal claches
- Only for visualization and lipid/protein ration and partition analyses. Not direcly adapted to subsequent molecular simulations


