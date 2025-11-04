# modisIceACI
Code for decomposing the shortwave effective radiative forcing from aerosol-cloud interactions (SW ERFaci) from ice clouds into components associated with the Twomey effect and IWP and ICF adjustments.

## References 
- Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: [Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part I: Cloud Radiative Kernels](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00248.1). J. Climate, 25, 3715-3735. 
    doi:10.1175/JCLI-D-11-00248.1.
- Zelinka, M. D., S. A. Klein, and D. L. Hartmann, 2012: [Computing and Partitioning Cloud Feedbacks Using 
    Cloud Property Histograms. Part II: Attribution to Changes in Cloud Amount, Altitude, and Optical Depth](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00249.1). 
    J. Climate, 25, 3736-3754. doi:10.1175/JCLI-D-11-00249.1.
- Zelinka, M.D., S.A. Klein, K.E. Taylor, T. Andrews, M.J. Webb, J.M. Gregory, and P.M. Forster, 2013: 
    [Contributions of Different Cloud Types to Feedbacks and Rapid Adjustments in CMIP5](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00555.1). 
    J. Climate, 26, 5007-5027. doi:10.1175/JCLI-D-12-00555.1.
 
## Input

The associated code requires the following inputs, including new effective radius x liquid water path cloud fraction joint histograms from the MODIS satellite simulator:

| Frequency | Name | Description | Unit | File Format |
|-----------|------|-------------|------|-------------|
| monthly mean | CLMODIS_IWPR | MODIS simulator ice-cloud fraction histograms | % | nc |
| monthly mean | FSDSC | Clearsky downwelling solar flux at surface | W/m^2 | nc |
| monthly mean | FSNSC | Clearsky net solar flux at surface | W/m^2 | nc |
| monthly mean | TS     | surface temperature | K     | nc            |
| monthly mean | SWkernel | SW cloud radiative kernel for ice clouds | W/m^2/% | nc |

The SW cloud radiative kernel is available to download at https://github.com/brandonduran/iceACI/tree/main/data

- ensmean_SW_kernel.nc: SW cloud radiative kernel developed using zonal mean temperature and humidity profiles averaged across control runs of five CMIP6-era climate models as input to the RRMTG radiation code. These are best for diagnosing feedbacks / forcing relative to a modeled pre-industrial climate state. Please refer to Wall et al. (2023) and Duran et al. (in prep) for details.

<!-- The inputs listed above are too large to host on Github. The data is available to download through the UCSD library digital collections: https://doi.org/10.6075/J0P26ZF1 -->

## Running the Notebook

Inside the data folder, two separate new folders should be created: "PI" and "PD". Each folder should contain four files, corresponding to the first four variables listed in the input table above for each simulation. The naming convention should be *VARNAME.nc*, eg *TS.nc*. 
<!-- These netCDF files are too big to host on Github, but are available through the UCSD library digital collections. -->
