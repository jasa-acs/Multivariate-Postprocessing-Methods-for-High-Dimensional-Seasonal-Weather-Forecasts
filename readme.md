# Multivariate Postprocessing Methods for High Dimensional Seasonal Weather Forecasts

# Author Contributions Checklist Form

## Data

### Abstract

Sea surface temperature observations (OISST dataset from NOAA), and raw
NorCPM forecast ensemble, for January/February 1985 – 2016, restricted
to longitudes between 90°W and 40°E, and latitudes between 0°N and 65°N.

### Availability

Restricting the dataset to two months and to a region containing only
8500 grid points, reduces the size of the raw dataset from roughly 4.5
GB to 71 MB. The data derived and saved in the course of postprocessing
has a size of roughly 14GB for the reduced dataset, and requires more
than 500 GB when the entire globe and all months are used.

Restricting the data spatially and temporally, while changing the values
of the derived scores, does not affect any of the conclusions derived in
Section 4 and 5 of the paper.

### Description

Code and data are available in form of an R package at

https://github.com/ClaudioHeinrich/pp.sst

This package contains the above described example data.

The full OISST observational dataset is publicly available at
ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc

In order to download the unprocessed NorCPM predictions, access to the
NIRD Norwegian data storage is required, see details at
<https://wiki.uib.no/norcpm>. The path to the run used in the paper is
projects/NS9039K/shared/norcpm/cases/NorCPM/NorCPM_V1c/ana_19800115_me_hindcasts.

Access to this run can be granted upon request, please contact the
authors of the article.

Both forecasts and observations are originally NetCDF files, but are in
a preliminary step coerced into an R data table. Both metadata and
instructions how to derive the data table from the original data are
contained in the README-file of the package and the package
documentation.

## Code

### Abstract

R code for all postprocessing methods described in the paper.

### Description

The R package pp.sst available at
*<https://github.com/ClaudioHeinrich/pp.sst>* contains both the code and
the data described above.

The code is under MIT license.

## Instructions for Use

### Reproducibility

The code reproduces all scores, PITs and rank histograms presented in
the paper (for the reduced dataset). The workflow follows along several
scripts available at
*<https://github.com/ClaudioHeinrich/pp.sst/tree/master/scripts>*

Detailed instructions for the reproduction are contained in the file
README-file available at *<https://github.com/ClaudioHeinrich/pp.sst>.*

The code was developed and tested on a server with 480 GB RAM and 30 CPU
cores, where one full run (all master scripts, considering the entire
provided test data set) took approximately 1 hour and 20 minutes to
complete. The memory usage never went above 2 %, and only two cores were
in fact used. The code will therefore also run on any standard laptop.
For shorter runtimes we recommend restricting the data to a smaller
spatial window. When run on Windows, parallelization needs to be turned
off. Instructions how to do this are contained in the README-file.

The reproduction relies on the data (both observations and predictions)
being already stored in a data table of the same format as the example
data table DT. The file scripts/various/data_processing.R shows how
such a data table can be constructed from NetCDF files.
