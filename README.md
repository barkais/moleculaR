
<img src="moleculaR_logo.png" width="600" height="120">

**An R package for the physical-organic chemist - designed for easy, chemical-intuition based molecular features extraction and statistical modeling.**
>**See moleculaR's [web page](https://barkais.github.io/) for detailed guides and  documenntation**
## Installation 

Currently, as the package is still under development (yet to be submitted to CRAN), installations proceed with the help of the `remotes` package.

It is highly recommended to use [RStudio]('https://posit.co/'), although applicable otherwise as well. Mac and Windows installations should not expect any weird behaviors with molecular visualizations and model plots, while the same can not be stated with linux installations (requires OpenGL). 

### Installation from Github 

First, install the `remotes` package from CRAN.
The `repos = getCRANmirrors()[1,"URL"]` addition helps installation in linux interactive sessions.

```
install.packages('remotes', repos = getCRANmirrors()[1,"URL"])
```

Once `remotes` is properly installed, use the `install_github` function to install `moleculaR`.
For convenience, load the package.

```
# Install
remotes::install_github('https://github.com/barkais/moleculaR.git')

# Load
library('moleculaR')
```

## Usage

NOTE:

With the current version - it is only possible to analyze correctly indexed Gaussian .log files.
Make sure the common substructure is numbered correctly.

### Featurization

**See `Getting Started with Examples` in articles for a detailed guide**

#### Step 1 - Pull information from Gaussain log files

Make sure all log files you want to analyze are in one location. 

```
# Run extractoR() (while in the .log files' directory)
extractoR()
```
Expect the following message, indicating everything worked fine. 

`Done!`

**This action results in a set of .feather files, which are light weight files that hold all the information `moleuclaR` needs. As we assume most users will prefer working on local machines, and to avoid transferring heavy log files to a local machine, it is best practice to install `moleculaR` on both the remote and the local, and to transfer only the resulting .feather files to the local.**

#### Step 2 - Unwrap feather files (either locally or remotely)

```
# Change working directory to folder containing the feather files resulted from running `extRactoR.auto()`
# Run unwRapper
unwRapper()
```

This command unwraps the feather files (compressed information from log files) and results in two new directories:

  1. Optimized_structures_xyz - containing xyz files extracted from the optimization.
  
  2. moleculaR_csv_files - A master directory that holds a sub directory for each of the molecules in the set. The sub directories contain csv files with all the information moleculaR needs, in human            readable format. 

Once finished, the function suggests to change your working directory. It is recommended to allow. 

#### Step 3 - Create moleculaR's input file

```
# Run moleuclaR input maker
moleculaR.Input.Maker()
```

**In case you opt to generate a 3D visualization, it is best practice to use one of the molecules found in the `Optimized_structures_xyz`.**

**Make sure to you do not save the input file in the moleculaR_csv_files folder. There should not be anything else in the moleculaR_csv_files folder.**

#### This step results in creating an input file in .txt format - containing the arguments for feature generation in `moleculaR()` - for future uses and documentation.

#### The resulting input file can now be passed to the main function `moleculaR()` as an argument. 

#### Step 4 - Get molecular features
```
# Use an input file generated with moleculaR.Input.Maker().
# Pass the input file location to moleculaR's main function `moleculaR()`.
# Indicate the path to the input file as an argument.
moleculaR('path/to/input/file.txt')
```

### Command Line Usage

It is possible to compute and extract all features directly from the command line in an interactive R session or as a part of a workflow. 

See `Command Line Usage` in articles for a detailed guide

### Modeling 

The package includes a set of linear and logistic regression modeling tools that allow users to input a dataset and get a cross validated (3-fold, 5-fold and LOO) model, along with nice plots. 

* To be able to use this functionality, user must add an _output_ or a _class_ column to the features.csv files resulted from running `moleculaR()`. The column *MUST* have the name 'output'/'class' (lower case intended). It is advised to make sure that the correct values are inserted in terms of order and naming - to prevent misleading model results.

See `Modeling` in articles for a detailed guide

### Cube files

In addition to the classic sterimol parameters, which are based on the xyz of an optimized structure and tabulated radii data, moleculaR also includes a cube file based analysis, along with an informative 3D visualization.

see `Cube functions` in articles for a detailed guide
