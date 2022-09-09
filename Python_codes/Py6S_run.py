# -*- coding: utf-8 -*-
"""
Authors: Rejane Paulino, Felipe Nincao, Thainara Munhoz
National Institute for Space Research (INPE) - Instrumentation Laboratory for Aquatic Systems (LabISA)

Last update: 08 september 2022
"""
# Libraries:

from Py6S import *
import pandas as pd
import numpy as np
import csv

#%%
# Input data:
    #aod550: Aerosol concentration at 550 nm
    #solar_Az: Solar azimuth angle (degree)
    #solar_Zn: Solar zenital angle (degree)
    #view_Zn: view zenith angle (degree)
    #alt_alvo: target altitude (km)

aod550 = 0.0781017
solar_Az = 43.9360618
solar_Zn = 37.66965103
view_Zn = 0.79983764
alt_alvo = 0.38607053

s = SixS()

# Atmospheric profile------------------------------------------------------------------------------------------:
s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

#Aerosol-------------------------------------------------------------------------------------------------------:
s.aot550 = aod550

#Aerosol profile-----------------------------------------------------------------------------------------------:
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)

# Ground Reflectance-------------------------------------------------------------------------------------------:
s.ground_reflectance = GroundReflectance.HomogeneousLambertian(np.array([[0.419, 0.427, 0.434, 0.442, 0.449, 0.456, 0.464, 0.471, 0.478,0.485, 0.493, 0.5  , 0.508, 0.515, 0.523, 0.531, 0.538, 0.546,0.555, 0.563, 0.571, 0.579, 0.588, 0.596, 0.605, 0.614, 0.623,0.632, 0.641, 0.651, 0.66 , 0.67 , 0.679, 0.689, 0.699, 0.709,0.719, 0.729, 0.739, 0.75 , 0.76 , 0.771, 0.781, 0.791, 0.802,0.813, 0.823, 0.834, 0.844, 0.855, 0.866, 0.877], [0.01445915, 0.01422249, 0.01406214, 0.01449518, 0.01592469,0.01761156, 0.01887088, 0.01983343, 0.02057033, 0.02136334, 0.02239862, 0.02424129, 0.02697532, 0.03056817, 0.03507325,0.04024491, 0.04571142, 0.05035407, 0.05309473, 0.0534503 ,0.05134707, 0.04744315, 0.04276923, 0.03773851, 0.03321408, 0.03015397, 0.02826281, 0.02778316, 0.02844564, 0.02861082,0.02545578, 0.02104192, 0.02031077, 0.0262125 , 0.03404532,0.03589206, 0.03120186, 0.023576  , 0.01819808, 0.01653794, 0.01646553, 0.01618563, 0.01658291, 0.01753715, 0.01833317,0.01814674, 0.01594022, 0.01273109, 0.01097039, 0.01016107,0.00947185, 0.00875745]]))

# Geometries of view and illumination--------------------------------------------------------------------------:
s.geometry = Geometry.User()
s.geometry.day = 4
s.geometry.month = 9
s.geometry.solar_z = float(solar_Zn)
s.geometry.solar_a = float(solar_Az)
s.geometry.view_z = float(view_Zn)

# Altitudes---------------------------------------------------------------------------------------------------:
s.altitudes = Altitudes()
s.altitudes.set_sensor_satellite_level()  # Set the sensor altitude to be satellite level.
s.altitudes.set_target_custom_altitude(alt_alvo)  # The altitude of the target (in km).

#%%
# Interpolation function -> 2.5 nm
def _interpolate_(dataframe_):
    import scipy.interpolate
    min_ = int(dataframe_.columns[0])+1
    max_ = int(dataframe_.columns[len(dataframe_.columns)-1])
    a = np.arange(min_,max_,2.5)
    new_srf = [] # Recebe os novos valores de SRF para cada banda
    for k in range(0,len(dataframe_)):
        y = dataframe_.iloc[k,:]
        x = dataframe_.columns
        y_interp = scipy.interpolate.interp1d(x, y, fill_value="extrapolate")
        new_srf.append(y_interp(a)) 
    arr = np.asarray(new_srf)
    new_bands = pd.DataFrame(arr, index = dataframe_.index, columns = a)
    return(new_bands)

#%%
# PRISMA - Spectral Response Function
frs = pd.read_excel('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA Spectral SWIR.xlsx',index_col=0)

frs_nova = _interpolate_(frs)
frs_nova_transposed = frs_nova.transpose()
frs = frs.transpose()

#%%
# Rodar o 6S para cada banda espectral    
dicionario = {}
waves_6s = []
waves_max_int = []

min_ = int(frs_nova_transposed.index[0])
max_ = int(frs_nova_transposed.index[len(frs_nova_transposed.columns)-1])

for a in range(0,frs_nova_transposed.shape[1]):
    
    max_value_index = frs_nova_transposed[frs_nova_transposed.iloc[:,a]==frs_nova_transposed.iloc[:,a].max()].index.values[0]
    waves_max_int.append(max_value_index)
    
    wv_value = frs_nova_transposed.columns.values[a]
        
    #srf_values = list(frs_nova_transposed.loc[379:1046.5,wv_value])
    srf_values = list(frs_nova_transposed.loc[min_:max_,wv_value])
    
    s.wavelength = Wavelength((min_/1000),(max_/1000),srf_values)
    
    s.run()
    
    print(s.wavelength)
    
    s.outputs.values['spherical_albedo'] = s.outputs.spherical_albedo.total
    s.outputs.values['co_transmittance_total'] = s.outputs.transmittance_co.total
    s.outputs.values['co2_transmittance_total'] = s.outputs.transmittance_co2.total
    s.outputs.values['oxyg_transmittance_total'] = s.outputs.transmittance_oxygen.total
    s.outputs.values['no2_transmittance_total'] = s.outputs.transmittance_no2.total
    s.outputs.values['ch4_transmittance_total'] = s.outputs.transmittance_ch4.total
    s.outputs.values['ozone_transmittance_total'] = s.outputs.transmittance_ozone.total
    s.outputs.values['water_transmittance_total'] = s.outputs.transmittance_water.total
    s.outputs.values['total_scattering_transmittance_upward'] = s.outputs.transmittance_total_scattering.upward
    s.outputs.values['total_scattering_transmittance_downward'] = s.outputs.transmittance_total_scattering.downward
    
    dicionario[str(int(max_value_temp))] = s.outputs.values

    
#%%
# Save output parameters in a json file
import json
a_file = open("G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA_04-09-21_resultados/6S_parameters_SWIR.json", "w")
json.dump(dicionario, a_file)
