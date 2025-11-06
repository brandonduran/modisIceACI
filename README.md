# modisIceACI

[![DOI](https://zenodo.org/badge/1086077645.svg)](https://doi.org/10.5281/zenodo.17544907)

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
  - Wall, C. J., Storelvmo, T., A. Possner, 2023: [Global observations of aerosol indirect effects from marine liquid clouds](https://acp.copernicus.org/articles/23/13125/2023/). Atmospheric Chemistry and Physics, 23, 13 125–13 141, doi:10.5194/acp-23-13125-2023.
  - Duran, B.M., Wall, C.J., Lutsko, N.J., Michibata, T., Ma, P.L., Qin, Y., Duffy, M.L., Medeiros, B. and Debolskiy, M., 2025: [A new method for diagnosing effective radiative forcing from aerosol–cloud interactions in climate models](https://acp.copernicus.org/articles/25/2123/2025/). Atmospheric Chemistry and Physics, 25(4), 2123-2146. doi:10.5194/acp-25-2123-2025.
 
## Input

The associated code requires the following inputs, including new effective radius x liquid water path cloud fraction joint histograms from the MODIS satellite simulator:

| Frequency | Name | Description | Unit | File Format |
|-----------|------|-------------|------|-------------|
| monthly mean | CLMODIS_IWPR | MODIS simulator ice-cloud fraction histograms | % | nc |
| monthly mean | FSDSC | Clearsky downwelling solar flux at surface | W/m^2 | nc |
| monthly mean | FSNSC | Clearsky net solar flux at surface | W/m^2 | nc |
| monthly mean | TS     | surface temperature | K     | nc            |
| monthly mean | SWkernel | SW cloud radiative kernel for ice clouds | W/m^2/% | nc |

Surface temperature is needed only for calculation of cloud feedbacks.

The SW cloud radiative kernel is available to download at https://github.com/brandonduran/iceACI/tree/main/data

- SWkernel_CTP250_iceflag3.nc: SW cloud radiative kernel developed using the zonal mean temperature and humidity profile from a CAM6 control run. These are best for diagnosing feedbacks / forcing relative to a modeled pre-industrial climate state. Please refer to Wall et al. (2023) and Duran et al. (2025) for details.

<!-- The inputs listed above are too large to host on Github. The data is available to download through the UCSD library digital collections: https://doi.org/10.6075/J0P26ZF1 -->

## Running the Notebook

Inside the data folder, two separate new folders should be created: "PI" and "PD". Each folder should contain four files, corresponding to the first four variables listed in the input table above for each simulation. The naming convention should be *VARNAME.nc*, eg *TS.nc*. Example data files are too big to host here, but will be made available through the UCSD library digital collections when the publication is accepted. Alternatively, all PPE output is accessible through the NCAR filesystem, or reach out to me for access.
