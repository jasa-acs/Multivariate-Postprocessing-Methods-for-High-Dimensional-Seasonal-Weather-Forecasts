
#################################################################
 Code repository accompanying the paper 'Multivariate postprocessing methods for high-dimensional seasonal weather forecasts' 
#################################################################

This package contains codes for the post-processing methods developped in the paper 'Multivariate postprocessing methods for high-dimensional seasonal weather forecasts'.
It requires R, version 3.4.1 or newer. Roughly 10GB RAM are minimal requirement for running the code on the full test data set. If less is available, reduce the size of the spatial window, see below.

The package can be installed by devtools::install_github('ClaudioHeinrich/pp.sst'). However, this command is included in the first script 01.master.setup.R

The package contains a data example of sea surface temperature forecasts issued by the Norwegian climate prediction model, as well as OISST observations, from 1985 to 2016, 
covering large parts of the Atlantic ocean (~70 MB). 
The OISST observational dataset is freely available at ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc
The unprocessed predictions by NorCPM require access to NIRD Norwegian, see details at https://wiki.uib.no/norcpm.
The link to the run used in the paper is /projects/NS9039K/shared/norcpm/cases/NorCPM/NorCPM_V1c/ana_19800115_me_hindcasts/.
If you are interested in getting access to the entire data, please contact one of the authors.
For more details on the data we refer to the main body of the paper as well as the data documentation in the R package.

The name of the example data table is DT, '?DT' provides further information. An important function to visualize the data is the function plot_diagnostic, see documentation.

The code was tested on a Linux x86 (64-bit) server with two Inter Xeon E5-2630 v3 processors, 16 cores, 32 threads and 256 GB memory, where one full run (all master scripts, considering the entire provided test data set) took approximately 1 hour and 20 minutes to complete. The memory usage never went above 2 %, and only two cores were in fact used. The code will therefore also run on any standard laptop. For shorter runtimes we recommend restricting the data to a smaller spatial window. This is done by changing the variables lat_box and lon_box in script 01.master.setup.R. When the code is run on Windows, parallelization needs to be turned off, which is done by setting the variable mc_cores to 1 in 01.master.setup.R.


######################
### Code structure ###
######################

The reproduction of the results follows along several scripts (in the scripts folder, i.e. at https://github.com/ClaudioHeinrich/pp.sst/tree/master/scripts). 
These scripts are divided into main scripts (named <number>.master.<description>.R) and side scripts
(named <masternumber>.side.<sidenumber>.<description>.R).
The numbering indicates dependencies, as described in the following. 

The main reason for defining this script structure is the complexity of the Post-processing methods, which can result in long computation times when large data sets are considered.
If everything would be merged into one massive script, this might result in unnecessary long computation times, since more often than not one is not interested in all parts of the analysis.
Splitting the scripts therefore makes it easier and faster to reproduce specific results by avoiding all unnecessary steps.
As argued in the article, computation time is usually not a major issue when dealing with seasonal forecasts. 
When you predict a month ahead, it is quite alright spending several hours and possible days to improve your forecast, unlike in short-range weather prediction.

The reproduction relies on the data (both observations and predictions) being already stored in a data table of the same format as the example data table DT. The file scripts/various/data_processing.R shows how such a data table can be constructed from netcdf files. 

Additionally, the script 'generate_plots' reproduces all plots shown in the paper. It essentially requires all other scripts to be run previously, see the script for more details.

###########################################################################################
##### The script 01.master.setup.R, handling multiple runs and the variable name_abbr #####
###########################################################################################

This script 01.master.setup.R is the first one that should be run. It installs the pp.sst package (which contains the example data table DT required for reproduction), and specifies parameters. 

Sometimes you might want to try different parameter settings (e.g. a different validation period) without overwriting the derived data from previous runs.
This is handled by the variable name_abbr, the abbreviated name of the run.
As an example, say you start a full run including global SST data and you set name_abbr to 'Full'. 
The script lets you set a directory name for saving data associated with this run. This name should (and does by default) contain the name_abbr. 
So you could set your save_dir for example to 'home/postprocessing/Full/'.
At some point halfway through the analysis you then realize that you're mostly interested in the SST of the Mediterranean sea. 
Moreover, you start to believe that the chosen training period is a bit too short. 
So you decide to start a new run with a shorter validation period (i.e. longer training period) and restrict locations to the Mediterranean sea, making the run much faster.
In order to not overwrite your previous 'Full' run you choose a new name_abbr, say 'Med_short_val'. 
Rerunning the 01.master.setup.R script with this name_abbr will now create a new directory for saving data, namely 'home/postprocessing/Med_short_val/'. 
Your previous 'Full' run will not be overwritten and you can return to it at any time.
When you use various runs with various specifications it's advisable to store a short description of the run in the 'description' variable.

The script then sets several variables and creates directories for data and plots, and does some small data modifications.
Finally, it saves its image as the file 'setup.R' in the specified save_dir. 
The 'setup.R' file is overwritten by each master script, see below, and contains the progress of your current run, specified by name_dir and the corresponding save_dir.
In our example above, in order to later return to the 'Full' run, you would therefore simply set name_abbr to 'Full', save_dir to 'home/postprocessing/Full/' and load 
the 'setup.R' file in save_dir, and you're good to go.

#####################################
##### All other master scripts ######
#####################################

All other master scripts require you first to specify your current run by setting name_abbr and specifying save_dir. 
save_dir should therefore depend on name_abbr such that it automatically updates when you change name_dir, and you don't need to change it manually every time.
In the example from the last section, you should therefore set your save_dir to paste0('home/postprocessing/',name_abbr,'/').

The script then loads the 'setup.RData' file in the specified save_dir. Thereafter it does its own specific piece of analysis, 
e.g. the second master script estimates the bias of the SST forecast as described in the paper.
At the end, it overwrites the 'setup.RData' file with its image.

The master scripts are designed to be run in succession. Jumping from the 1st to the 3rd master script, for example, will not work, 
since required variables have not yet been defined, datasets have not yet been created, etc. The progress of the current run is stored in the script_counter variable,
which is particular helpful when you switch between runs, as you don't need to remember where you left the old run.

############################
##### The side scripts #####
############################

The side scripts contain parts of the analysis that are not fundamental and may be skipped. Like the master scripts, they first specify name_abbr and load the corresponding setup.RData file. 
They have a masternumber, which indicates the master script they depend on, and a side number which is merely a suggested ordering. 
For example, the side script 03.side.06.PIT.plots.R requires you to have previously run the third master script, but you can run it before 03.side.03.permutation.tests.univariate.R.
Moreover, you can run the 4th master script without running either one of them, and you can still run both of them later, say after you finish the 5th master script.
 
Both master and side scripts have a header that clears the workspace, loads packages, sets name_abbr, and loads the setup.RData file. 
If you run several scripts in succesion you may want to skip this, since the setup file can grow quite big if e.g. the entire globe is considered.




