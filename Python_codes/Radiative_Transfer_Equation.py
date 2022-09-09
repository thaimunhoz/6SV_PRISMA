# -*- coding: utf-8 -*-
"""
Authors: Rejane Paulino, Felipe Nincao, Thainara Munhoz
National Institute for Space Research (INPE) - Instrumentation Laboratory for Aquatic Systems (LabISA)

Last update: 08 september 2022
"""

# Library:
import numpy as np
from osgeo import gdal
import numpy as np
import pandas as pd
import json
import os
from PIL import Image

#%%
# Excel file containing the central wavelenght, fwhm, name of each band for the specific image
# IMPORTANT -> The number of bands must be equal to the spectral response function
bands_prisma = pd.read_excel(
    'G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/prisma_bands_data.xlsx', header=0, index_col=0)
bands_prisma_acolite = bands_prisma.iloc[:, 1]

# TOA Reflectance -> each band must be saved in a unique folder
path_to_toa_reflectance = 'G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA_04-09_rtoa/rhot_'

# Read the JSON file for the VNIR
with open('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA_04-09-21_resultados/6S_parameters_VNIR.json') as f:
    data = f.read()

dicionario_VNIR = json.loads(data)

# Read the JSON file for the SWIR
with open('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA_04-09-21_resultados/6S_parameters_SWIR.json') as f:
    data = f.read()

dicionario_SWIR = json.loads(data)

dicionario_completo = {**dicionario_VNIR, **dicionario_SWIR}

dic_aux = dicionario_completo

a = 0
for k in list(dic_aux.keys()):
    dicionario_completo[str(bands_prisma_acolite[a])] = dicionario_completo.pop(str(k))
    a = a + 1

#%%
# Apply the radiative transfer equation in each RTOA image considering the parameters from the 6S
for a in bands_prisma_acolite.values:

    string_band = path_to_toa_reflectance + str(a) + ".tiff"
    band_prisma = gdal.Open(string_band, gdal.GA_ReadOnly)
    array_prisma = band_prisma.GetRasterBand(1).ReadAsArray()

    # Atmospheric modeling:
    atmospheric_modeling_py6S_Band = dicionario_completo[str(a)]

    # Atmospheric Correction -> Equation Vermote et al. (2016):
    tg_OG_co = float(atmospheric_modeling_py6S_Band['co_transmittance_total'])
    tg_OG_c02 = float(
        atmospheric_modeling_py6S_Band['co2_transmittance_total'])
    tg_OG_o2 = float(
        atmospheric_modeling_py6S_Band['oxyg_transmittance_total'])
    tg_OG_no2 = float(
        atmospheric_modeling_py6S_Band['no2_transmittance_total'])
    tg_OG_ch4 = float(
        atmospheric_modeling_py6S_Band['ch4_transmittance_total'])

    # Total transmission of Other Gases -> Tg_OG:
    Tg_OG = float(tg_OG_co * tg_OG_c02 * tg_OG_o2 * tg_OG_no2 * tg_OG_ch4)

    # Total transmission of the Ozone -> Tg_O3:
    Tg_O3 = float(atmospheric_modeling_py6S_Band['ozone_transmittance_total'])

    # Total transmission of the Water Vapor -> Tg_H2O:
    Tg_H20 = float(atmospheric_modeling_py6S_Band['water_transmittance_total'])

    # Total transmittance upward (Rayleigh + Aerosol) -> T_upward:
    T_upward = float(
        atmospheric_modeling_py6S_Band['total_scattering_transmittance_upward'])

    # Total transmittance downward (Rayleigh + Aerosol) -> T_downward:
    T_downward = float(
        atmospheric_modeling_py6S_Band['total_scattering_transmittance_downward'])

    # Total transmission of the atmosphere -> T_atm:
    T_atm = float(T_upward * T_downward)

    # Atmosphere intrinsic reflectance -> p_atm:
    p_atm = float(
        atmospheric_modeling_py6S_Band['atmospheric_intrinsic_reflectance'])

    # Atmosphere spherical albedo -> s_atm:
    s_atm = float(atmospheric_modeling_py6S_Band['spherical_albedo'])

    # Equation (Ref: Vermonte et al., 1997)
    f1 = Tg_OG * Tg_O3 * Tg_H20 * T_atm
    f2 = p_atm / (Tg_H20 * T_atm)
    parte_1 = np.divide(array_prisma, f1)
    parte_2 = np.subtract(parte_1, f2)
    parte_3 = np.multiply(parte_2, s_atm)
    parte_4 = np.add(parte_3, 1)
    p_s = parte_2 / parte_4

    rrs = np.divide(p_s, np.pi) # Calculate de Rrs

    rrs[rrs == rrs[0][0]] = 0

    img_tiff = Image.fromarray(rrs,mode='F')
    # Save each surface reflectance image in a separate file
    name = 'C:/Users/thain/Documents/6S_run/rhos_' + str(a) + '.tiff'
    img_tiff.save(name,'TIFF')
    
    band_prisma = None

    print(a)

print('FIM')

# %%
# Join all the bands in a stack file:
file_list = ['C:/Users/thain/Documents/6S_run/rhos_' +
             str(k) + '.tiff' for k in bands_prisma_acolite.iloc[0:63]] #Considering only the VNIR (bands 0-63)

files_string = " ".join(file_list)

command_to_join = "gdalbuildvrt -separate stack_6s.vrt " + files_string

os.system(command_to_join)

os.system('gdal_translate stack_6s.vrt PRISMA_6s_rrs_VNIR.tif')

#%%
#Clean the folder
for k in bands_prisma_acolite.values:
    if os.path.exists('C:/Users/thain/Documents/6S_run/rhos_' + str(k) + '.tiff'):
        os.remove('C:/Users/thain/Documents/6S_run/rhos_' + str(k) + '.tiff')