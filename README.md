---
output: html_document
---

<img src="man/figures/moleculaR_logo.png" width="600" height="120">

**An R package for the physical-organic chemist - designed for easy, chemical-intuition based molecular features extraction and statistical modelling.**

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
Make sure the common substructure is numbered correctly, **and that there is a copy of the simplest molecule (i.e. the molecule with the smallest substituent) named basic.log**

### Featurization

**See `Getting Started with Examples` in articles for a detailed guide**

#### Step 1 - Pull information from Gaussain log files

Make sure all the log files you want to analyze are in one location. 

```
# Run extRactoR (no need to change working directory - extRactoR does it for you)
extRactoR()
```
Running `extRactoR()` initiates a user interface, which allows for the selection of the log files directory, 
naming the resulting directory (created by the function) and sets the working directory to that. 

Expect the following message, indicating everything worked fine. 

`Done!`
` `      
`Current working directory is set to path/to/example_log_files/Name_Chosen_By_User`

**This action results in a set of .feather files, which are light weight files that hold all the information `moleuclaR` needs.** **As it is reasonable to assume that most users will prefer working on local machines, it is still recommended** **avoiding thransferring heavy log files to locals, and it is best practice to install `moleculaR` on both the remote and the** **local, and to transfer only the resulting .feather files to the local. **

#### Step 2 - Unwrap feather files (either locally or remotely)

```
# Run unwRapper (no need to change working directory - unwRapper does it for you)
unwRapper()
```

This command unwraps the feather files (compressed information from log files) and results in two new directories:

  1. Optimized_structures_xyz - containing xyz files extracted from optimization.
  
  2. moleculaR_csv_files - A master directory that holds a sub directory for each of the molecules in the     set. The sub directories contain csv files with all the information moleculaR needs, in human            readable format. 

Once finished, the function suggests to change your working directory. It is recommended to allow. 

#### Step 3 - Get molecular features

```
# Run moleuclaR (assuming you allowed unwRapper to change your working directory)
moleculaR()
```

**In case you opt to quit and generate a 3D visualization, it is best practice to use one of the molecules found in the** **`Optimized_structures_xyz`. After the visualization was generated, leave it open and run** 
**`moleculaR()` again.**

###### This step results in the creation of two files:
######  1. The output csv file - containing all wanted features.
######  2. An inputs file in .RData format - containing all user inputs to `moleculaR()` - for future uses and documentation.
#####  3. In cases where the exact same instructions should be passed to `moleculaR()`, you can use `moleculaR.input()` instead. 

```
# To use an input file resulting from moleculaR(), do one of the following:
#   1. Explicitly indicate the path to the input file as an argument
moleculaR.input('path/to/input/file.RData')
#   2. Run the function with no arguments and choose the file
moleculaR.input()
```

#### Command Line Usage

It is possible to compute and extract all features directly from the command line in an interactive R session or as a part of a workflow. 

See `Command Line Usage` in articles for a detailed guide

### Modeling 

The package includes a set of linear regression modeling tools that allow users to input a dataset and get a cross validated (3-fold, 5-fold and LOO) model, along with a nice plot. 

* To be able to use this functionality, user must add an _output_ column to the features.csv files resulted from running `moleculaR()` or `moleculaR.input()`. The column *MUST* have the name 'output' (lower case intended). It is highly advised to pay close attention that the correct experimental observations are attached, specifically - pay attention to molecule names.

See `Modeling` in articles for a detailed guide

### Cube files

In addition to the classic sterimol parameters, which are based on the xyz of an optimized structure and tabulated radii data, moleculaR also includes a cube file based analysis, along with an informative 3D visualization.

see `Cube functions` in articles for a detailed guide
