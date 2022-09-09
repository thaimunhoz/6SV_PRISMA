# 6SV_python in PRISMA images
Authors: Rejane Paulino, Felipe Nincao, Thainara Munhoz

National Institute for Space Research (INPE) - Instrumentation Laboratory for Aquatic Systems (LabISA)

Application of the 6S Radiative Transfer Model using the Py6S python library on PRISMA images. It allows you to run 6S simulations using Py6S and apply the simulation to the radiative transfer equation for atmosphere correction on the PRISMA L1 images.

**PY6S_run script**

Input data:
  - Aerosol concentration at 550nm
  - Atmospheric profile (**AtmosProfile**) → It considers the characteristics associated with this profile such as the amount of water vapour, pressure and ozone and incorporates it into the model. You can use a ready-made profile that is chosen based on latitude (`AtmosProfile.Tropical, AtmosProfile.MidlatitudeSummer, AtmosProfile.MidlatitudeWinter`), you can enter specific values ​​for water vapor and ozone*,* or you can to define a profile from radiosonde measurements.
  - Aerosol model (**AeroProfile**) → Defines the distribution, particle size, physical characteristics, etc. It can be considered a predefined type or set specific information.
  - Surface reflectance (**GroundReflectance**) → Produces strings for the input file for a number of different ground reflectance scenarios (Lambertian Homogeneous or BRDF, or Lambertian Heterogeneous). One can consider a constant spectral reflectance (a single value); a model reflectance spectrum for the study environment; or you can use some surfaces already inserted in the model (`GroundReflectance.GreenVegetation`).
  - Geometry information → Day, month, target altitude, solar zenith and azimuth angle, and zenith angle of sight.
  - Spectral response function of the sensor of interest in an Excel file format.
  
Output data:
  - Two JSON files (VNIR and SWIR) containing a dictionary with all the parameters necessary for the application o the Radiative Transfer Equation
  
**Radiative_Transfer_Equation scrypt**

Input data:
   - Prisma_band_data → Excel file containing: central wavelength of each band, fwhm, and the name of the band (this name must be associated with the name of the saved R_TOA images;
  - TOA Reflectance Images → One image for each band, saved separately in a single folder
  - JSON files containing Py6S VNIR and SWIR parameters
  
Output data:
  - Surface Remote Sensing Reflectance image.
