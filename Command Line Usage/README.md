# Command Line Usage - Guide

**This guide demonstrates the use of moleculaR's stand-alone functionalities.**

Note that there are some additional features that were not included in the `moleculaR()` function and are only available with interactive use. 

```r
# Load moleculaR
library(moleculaR)
```
***

## Molecule Visualization 

moleculaR's visual aid is accessible from the command line (for Mac and Windows installations, Linux users may be required to install openGL first).

```r
# To generate a plot of an xyz file, use:

plot_molecule(xyz_file)
```
***

>**All functions (but `steRimol.xyz` functions) require users to set the working directory to _moleculaR_csv_files_**

**Throughout the text, command-line usage will be regarded as CLU**

**Generally - moleculaR's CLU functions have suffixes .df and .multi:**
  
  **.df functions return data frames with wanted information**
  
  **.multi functions are vectorized versions of .df functions, applied to more than one input.**
  
**If there's any logic to a having a .multi function it will be there, but when there isn't, it won't.**
  
  

## steRimol
 
moleculaR's sterimol implementation for the calculation of sterimol values based on xyz structures. 

There are two distinct versions of steRimol, which depend on the task users wish to perform:
  
  1. `steRimol.df()` and `steRimol.multi()` - operate on moleculaR's csv files

This is the working function behind `moleculaR()` and the preferred function of use, as the rest of the package functions demand having the working directory set to _moleculaR_csv_files_.   

  2. `steRimol.xyz.df()` and `steRimol.xyz.multi()` - operate directly on xyz files 

Due to popular usage of sterimol values and the fact that the optimized structure is all it requires, this function allows users to have them calculated fast from any given xyz file. 

CLU for both cases is identical, and depends solely on the type of files you have in your current working directory.

#### TL;DR

```r
#Classic Use - for primary axis 1-2 use:

steRimol.multi('1 2') # for moleculaR's csv files

steRimol.xyz.multi('1 2') # for a set of xyz files

#For more than one primary axis, say 1-2 and 1-3, use:

steRimol.multi(c('1 2', '1 3')) # for moleculaR's csv files

steRimol.xyz.multi(c('1 2', '1 3')) # for a set of xyz files

#See documentation - in R console:
?steRimol.multi()
?sterimol.xyz.multi()
```

#### Function Structure

The structure of `steRimol.multi()` is presented, but is identical to that of `steRimol.xyz.multi()`
```r
steRimol.multi(coordinates_vector,
                CPK = T,
                only_sub = T,
                drop = NULL)
```

CPK, only_sub and drop all have default values, which allow immediate use with only primary axis supplied by users. 

##### Arguments:

  **coordinates_vector:** 
  
  An R vector of strings. Each string denotes a primary axis, with the structure 'atom1 atom2', where the 3D vector of atom1 to atom2 will define sterimol's primary axis. 
  
  The function expects at least one primary axis, when only one primary axis is wanted, user defines it with the argument:
  
  `'atom1 atom2'`
  
  To use multiple axes, user defines the set of axes by the argument:
  
  `c('atom1 atom2', 'atom3 atom4')`
  
  For example: 
  1. say a user wants to compute sterimol values for axis 1-2
  
  ```r
  # explicitly 
  steRimol.multi(coordinates_vector = '1 2')
  
  # implicitly 
  steRimol.multi('1 2')
  ```

    
  2. say a user wants to compute sterimol values for axes 1-2 and 1-3
  
  ```r
  # explicitly 
  steRimol.multi(coordinates_vector = c('1 2', '1 3'))
  
  # implicitly 
  steRimol.multi(c('1 2', '1 3'))
  ```

  **CPK:**
  
  Users can avoid default use of CPK radii and apply Pyykko's covalent radii values by setting `CPK = F`.
  
  
  For example:
  
  say a user wants to compute sterimol values for axes 1-2 and 1-3, with Pyykko's covalent radii values
  
  ```r
  # explicitly 
  steRimol.multi(coordinates_vector = c('1 2', '1 3'), CPK = F)
  
  # implicitly 
  steRimol.multi(c('1 2', '1 3'), F)
  ```
  
  **only_sub:**
  
  A default setting that makes sure that sterimol values are given only for the substituent in mind. steRimol uses a graph-based search for all atoms that stem from the primary axis, and uses the remaining substructure alone. 
  
  <center><img src="figures/sterimol_CLU.png" width="600" height="220"></center>
  
  
Users can avoid this step by setting `only_sub = F`

  ```r
  # CPK is either T or F
  
  # explicitly 
  steRimol.multi(coordinates_vector = c('1 2', '1 3'), CPK, only_sub = F)
  
  # implicitly 
  steRimol.multi(c('1 2', '1 3'), CPK, F)
  ```

  **drop:**
  
  An optional dropping of a substructure of the substituent in mind. Users can ask sterimol to drop all atoms that are directly connected from atom 'A' onward. 
  
  <center><img src="figures/sterimol_drop.png" width="500" height="250"></center>

In the example above, drop was set to 15. Note that for `drop`, user input is an integer and not a character. It is enough to indicate the first atom on the substructure that should be dropped. Multiple substructures can be dropped by using `drop = c(atom, another atom, ...)`.

  ```r
  # CPK is either T or F
  # only_sub is set to T, but can be F as well
  
  # explicitly 
  steRimol.multi(coordinates_vector = '3 1', CPK, only_sub = T, drop = 15)
  
  # implicitly 
  steRimol.multi('3 1', CPK, T, 15)
  ```
  
## Charges

Users can extract NBO, Hirshfeld, and CM5 charges, and differences between those extracted with the use of `nbo.df()`, `hirsh.df()`, and `cm5.df()`. The funcions work the same way as `nbo.df()` explained below. 

#### TL;DR

```r
# To extract NBO charges for atoms 1,2 and 3, use:

# explicitly
nbo.df(atom_indices = '1 2 3')

# implicitly 
nbo.df('1 2 3')

# and to add differences between those, say - 1 with 2, and 1 with 3, use:

# explicitly
nbo.df(atom_indices = '1 2 3', difference_indices = '1 2 1 3')

# implicitly 
nbo.df('1 2 3', '1 2 1 3')

# See documentation - in R console:
?nbo.df()
```

#### Function Structure

```r
nbo.df(atom_indices,
       difference_indices = NA)
```
The function has no default values, but differences is not used unless requested by user. 

##### Arguments:

  **atom_indices:**
  
  Indices of atoms of interest. Input is given as a string of indices, separated by a space. 
  
  For example, to get NBO charges of atoms 1,2,3 and 15, the argument should be `'1 2 3 15'`.
  
  **difference_indices:**
  
  Input is given as a string of indices (must be an even number of indices). The function interprets each pair (first & second, third & fourth, and so on) as a pair that should have its difference calculated. 
  
  For example, to get the differences between 1 and 2, and 1 and 3, the argument should be `'1 2 1 3'`
  



## Dipole Moments

See `Features - Definitions and Practice` for a detailed definition of the dipole moment as a feature.
 
 
### Gaussian's Dipole Moment

Users can extract dipole moment (components and magnitude) as it is given by Gaussian, or they may apply algebraic manipulations to it with chemical rationale in mind.

Change of coordinate system - the dipole's vector components with respect to a coordinate system of choice allows the use of the these components as chemical descriptors. [independent features](%22https://www.nature.com/articles/s41557-019-0258-1%22)

To change the coordinate system, users must define the new coordinate system with the following:

1)  A new origin - can be an atom or the centeroid of several atoms.
2)  A new y axis direction - an atom that defines the direction at which the y axis points.
3)  A new xy plane - an atom that defines the plane on which the x-axis lies. This atom does not determine the x-axis. Instead, it defines, together with the y-axis, the xy plane.

Users can extract the dipole moment for multiple coordinate systems.

#### TL;DR

```r        
# Dipole Moment = DM
# Coordinate System = CS

# To extract DM as it is given in Gaussian, use:

dip.gaussian.multi() # with no arguments

# To extract the DM, with a change of CS, such that atom 1 will be the origin,
# atom 2 will set the y direction, and atom 3 will set the xy plane, use:

# explicitly
dip.gaussian.multi(coor_atoms = '1 2 3')

# implicitly 
dip.gaussina.multi('1 2 3')

# To use the centeroid of a substructure (e.g. the center of atoms 3-8) 
# as the origin, with atom 1 defining the y direction and atom 4 defining the
# xy plane, use: 

# explicitly
dip.gaussian.multi(coor_atoms = '3 4 5 6 7 8 1 4')

# implicitly 
dip.gaussina.multi('3 4 5 6 7 8 1 4')

# To extract the DM using multiple CSs, input a vector of arguments.
# The '' input will provide the default Gaussian vector. 

# explicitly
dip.gaussian.multi(coor_atoms = c('3 1 4', '3 4 5 6 7 8 1 4', ''))

# implicitly 
dip.gaussian.multi(c('3 1 4', '3 4 5 6 7 8 1 4', ''))

# See documentation - in R console:
?dip.gaussian.df()
?dip.gaussian.multi()
```

#### Function Structure

```r         
dip.gaussian.multi(
  coor_atoms = c("")
)
```

##### Arguments:

**coor_atoms:**

Defaults to "", with no change of coordinates.

To change CS, user inputs indices of 3 or more atoms:

Atom indices, separated by a space - define:

i)  atom 1-3rd last - define the origin (the centroid of one or more atoms)
ii) atom 2nd last  - the y direction
iii) last atom - the xy plane
  
  
### Substructure isolated Dipole Moment (NPA Based)

In addition to Gaussian's dipole moment, users can extract the NPA (NBO) based dipole moment of a substructure. This option allows for a depiction of the charge distribution and direction of charge densities in a specific substructure. This classifies as **very** esoteric, and requires design in order to make sense. 

Note - the dipole moment calculated by this function takes ONLY the atoms requested into account, namely, it calculates the dipole moment created by the chosen atoms as if they were isolated from the molecule, and for that it uses the NPA charges of the atoms. 

For instance, say a user wishes to calculate the dipole moment created by the 6 carbons of a ring that is part of a molecule, as well as the ring including hydrogens and other substituents (that are consistent throughout the set):

<center><img src="figures/npa_dm.png" width="600" height="476"></center>

Users have to define the coordinate system as starting point. The coordinate system will not change, so it makes sense to use more than one substructure ONLY when the same coordinate system will be meaningful (as in the figure above). 

In this instance 

#### TL;DR

```r
# To extract the dipole moments as they appear in the figure above,
# the coordinate system should be defined as follows:
# 1. The origin is set to the moddle point between atoms 15 and 18
# 2. The y direction is set to atom 1 (or 15, the same)
# 3. The xy plane defining atom is set to 20 (or 19, the same)

# explicitly
npa.dipole.subunit.multi(coor_atoms = '15 18 1 20',
                         subunits_inputs_vector = c('15 16 17 18 19 20', # for the ring on the right
                         '15 16 17 18 19 20 21 22 23 24 25 26 1')) # for the ring on the left

# implicitly
npa.dipole.subunit.multi('15 18 1 20',
                        c('15 16 17 18 19 20', # for the ring on the right
                        '15 16 17 18 19 20 21 22 23 24 25 26 1')) # for the ring on the left

```
#### Function Structure

```r
npa.dipole.subunit.multi(coor_atoms,
                         subunits_inputs_vector)
```

##### Arguments:

**coor_atoms:**

To change CS, user inputs indices of 3 or 4 atoms:
  
  1. 3 atoms, separated by a comma - define: 
          
        i) atom 1 - the origin 
        ii) atom 2 - the y direction
        iii) atom 3 - the xy plane
        
  For example the input `'1 2 3'` will transform the coordinate system such that 1 is the origin, 2 is the y direction and 3 defines the xy plane.
          
  2. 4 atoms, separated by a comma - define:
  
        i) atoms 1 and 2  - center point between the two is the origin 
        ii) atom 3 - the y direction
        iii) atom 4 - the xy plane
        
  For example the input `'1 2 3 4'` will transform the coordinate system such that the center point between 1 and 2 is the origin, 3 is the y direction and 4 defines the xy plane.
  
**subunits_inputs_vector:**

vector of atoms characters. 
  
  For example: `c('10 11 12 13', '4 5 6 7 8 9')`
  
  Users have to define the coordinate system as starting point. The coordinate system will not change, so it makes sense to use more than one substructure ONLY when the same coordinate system will be meaningful.
  
  
## Vibrational Frequencies 

See `Features - Definitions and Practice` for a detailed definition of vibrational features and the algorithmic reasoning behind their extraction. 

### Bond Stretch

Users can directly extract the stretching vibrational frequencies of bonded pairs of atoms. Note that atoms must be bonded, and that submitting non-bonded atoms will crash the program with an error message.

#### TL;DR

```r
# To extract vibrational frequencies of bond 1-27 and 11-12, use:

# explicitly 
stretch.vib.df(atom.pairs = '1 27 11 12')

# implicitly
stretch.vib.df('1 27 11 12')
```

#### Function Structure

```r
stretch.vib.df(atom.pairs,
               threshold = 1350)
```

threshold is defined as the minimal wave number allowed to be assigned to a certain vibrational mode, as below this value movements are not distinct enough (fingerprint region).
Users may apply different thresholds, but should expect weird behavior with any value lower than 1350 cm^-1^

##### Arguments:

**atom.pairs:**

Input is given as a string of indices (must be an even number of indices). The function interprets each pair (first & second, third & fourth, and so on) as a pair that should have its bond's stretching frequency extracted. 

For example, to extract the frequencies of bonds 1-27 and 11-12, use: `'1 27 11 12'`.

**threshold:**

minimal wave number value (above fingerprint region). A numeric value. 



### Ring Vibrations

Ring vibrations are defined with one of the ring's six atoms, that defines the "Primary" atom of the ring.

<center><img src="figures/rings.png" width="600" height="382"/></center>

Users arbitrarily choose a "primary" atom - it is most convenient to choose the first atom on the ring, that connects the ring to the common substructure, though it really doesn't matter.

Note that it is possible to extract vibrations for as many rings wanted.

In this instance, our input would be `c(3, 15)` with `3` being the primary for the ring on the right, and `15` for the ring on the left.

#### TL;DR

```r       
# To extract ring vibrations for two rings, indexed 3-8 and 15-20
# (having the structures in the figure above as reference). Define, for instance,
# atoms 3 and 15 as primaries -> get the following inputs:

# Ring 1: primary = 3
# Ring 2: primary = 15

# The command will be:

# explicitly 
ring.vib.multi(inputs_vector = c(3, 15))

# implicitly
ring.vib.multi(c(3, 15))
```

#### Function Structure

```r         
ring.vib.multi(inputs_vector)
```

##### Arguments:

**inputs_vector:**

vector of characters of single atom indices.


### Bending Vibrations

Users can directly extract the bending vibrational frequencies of atoms that share a center (a third atom bonded to both of them). Note that atoms must have this center atom. This will work with any two atoms that stand this condition, namely, even if they do not actually move in what we would define as a bending motion, so we advise on using this with attention to details, and in proper settings which are chemically intuitive.  

#### TL;DR

```r
# This is a dummy example that is not relevant to the example files supplied, as they do not have 
# any relevant vibrational modes. 

# To extract vibrational bending frequencies of atoms 1 and 3, both bonded to the center atom 2,
# and similarly for atoms 4 and 6, that share the center 5, use:

# explicitly 
bend.vib.df(atom_pairs = '1 3 4 6')

# implicitly
bend.vib.df('1 3 4 6')
```

#### Function Structure

```r
bend.vib.df(atom_pairs)
```

##### Arguments:

**atom_pairs:**

Input is given as a string of indices (must be an even number of indices). The function interprets each pair (first & second, third & fourth, and so on) as a pair that should have its bending frequency extracted. 

## Angles

Users can extract angles between 3 atoms and dihedrals created with 4 atoms. The use of 3 atoms (a1, a2, a3) essentially measures the angle created by the bonds a1-a2 and a2-a3. The use of 4 atoms (a1, a2, a3, a4) measures the angle created by the planes of a1-a2-a3 and a2-a3-a4. 

Users can extract how many angles and dihedrals as wanted. Note that there are no limitations on the logic behind the calculation, namely, angles and dihedrals that actually mean something - any series of 3 or 4 atoms that exist in the molecule will return something, whether it makes chemical sense or not. 

#### TL;DR

```r
To extract the angle created by bonds 1-2 and 2-3, use:

# explicitly 
mol.angles.multi(atoms_vector = '1 2 3')

# implicitly
mol.angles.multi('1 2 3')

# To extract more than one value, use a vector of inputs. For example,
# to extract the angle created by bonds 1-2 and 2-3, and the dihedral of
# planes 4-3-1 and 3-2-1, use:

# explicitly 
mol.angles.multi(atoms_vector = c('1 2 3', '4 3 2 1'))

# implicitly
mol.angles.multi(c('1 2 3', '4 3 2 1'))
```

#### Function Structure

```r
mol.angles.multi(atoms_vector)
```

##### Arguments:

**atoms_vector:**

A vector of strings of either 3 or 4 atom indices. Users can mix between kinds and the result will reflect the choice. 

## Distances

Users can extract distances between pairs of atoms, without any limitation. 

#### TL;DR

```r
# To extract the distance between atoms 1-2 and 1-15, use:

# explicitly
atoms.distance.df(atom_pairs = '1 2 1 3')

# implicitly
atoms.distance.df('1 2 1 3')
```

#### Function Structure

```r
atoms.distance.df(atom_pairs)
```

##### Arguments:

**atom_pairs:**

Input is given as a string of indices (must be an even number of indices). The function interprets each pair (first & second, third & fourth, and so on) as a pair that should have its distance extracted. 

  
## Polarizabilty

Very simple usage:

```r
polar.df() # no arguments
```



