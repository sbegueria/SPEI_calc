# SPEI_calc

*A Standardized Precipitation-Evapotranspiration Index (SPEI) calculator written in C*

 
## Objective

The program calculates time series of the Standardised Precipitation-Evapotransporation Index (SPEI). 
 
## Use

The program is executed from the Windows console.
From an input data file containing monthly time series of precipitation and mean temperature, plus the geographic coordinates of the observatory, the program computes the SPEI accumulated at the time interval specified by the user, and generates a new data file with the SPEI time series.
It is easy to create a batch script for automating  the calculation of the SPEI over a large number of observatories or for several accumulated periods.
 
## IMPORTANT NOTICE

This software has been superseded by the SPEI package for R, which is currently the only development branch. Please, have a look at (CRAN)[http://cran.r-project.org/web/packages/SPEI/index.html], and (GitHub)[https://github.com/sbegueria/SPEI].
