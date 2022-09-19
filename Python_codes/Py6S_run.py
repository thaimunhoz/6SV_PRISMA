# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:19:12 2022

@author: thain
"""

from Py6S import *
import pandas as pd
import numpy as np
import csv

#%%
# Variáveis de entrada
aod550 = 0.060  # PRISMA 03/10/2021
solar_Az = 55.23819 # Angulo solar azimutal PRISMA 03/10/2021
solar_Zn = 27.972664 # Angulo solar zenital PRISMA 03/10/2021
view_Zn = 0.83701292 # Angulo azimutal do sensor PRISMA 03/10/2021

#aod550 = 0.0781017
#solar_Az = 43.9360618
#solar_Zn = 37.66965103
#view_Zn = 0.79983764

alt_alvo = 0.38607053 #km

#%%
# Função que permite o acesso ao modelo 6S
s = SixS()
#s.outputs.write_output_file('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/output_py6s.txt')

#%%
# Perfil atmosférico:
# Considera as características associadas a este perfil (vapor d'água e ozonio) e incorpora no modelo
s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

#Aerossol:
s.aot550 = aod550
#Modelo de aerossol a ser utilizado (define a distribuição, tamanho das partículas, caracteristicas físicas, etc)
s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)

# Ground Reflectance-------------------------------------------------------------------------------------------:
s.ground_reflectance = GroundReflectance.HomogeneousLambertian(np.array([[0.419, 0.427, 0.434, 0.442, 0.449, 0.456, 0.464, 0.471, 0.478,0.485, 0.493, 0.5  , 0.508, 0.515, 0.523, 0.531, 0.538, 0.546,0.555, 0.563, 0.571, 0.579, 0.588, 0.596, 0.605, 0.614, 0.623,0.632, 0.641, 0.651, 0.66 , 0.67 , 0.679, 0.689, 0.699, 0.709,0.719, 0.729, 0.739, 0.75 , 0.76 , 0.771, 0.781, 0.791, 0.802,0.813, 0.823, 0.834, 0.844, 0.855, 0.866, 0.877], [0.01445915, 0.01422249, 0.01406214, 0.01449518, 0.01592469,0.01761156, 0.01887088, 0.01983343, 0.02057033, 0.02136334, 0.02239862, 0.02424129, 0.02697532, 0.03056817, 0.03507325,0.04024491, 0.04571142, 0.05035407, 0.05309473, 0.0534503 ,0.05134707, 0.04744315, 0.04276923, 0.03773851, 0.03321408, 0.03015397, 0.02826281, 0.02778316, 0.02844564, 0.02861082,0.02545578, 0.02104192, 0.02031077, 0.0262125 , 0.03404532,0.03589206, 0.03120186, 0.023576  , 0.01819808, 0.01653794, 0.01646553, 0.01618563, 0.01658291, 0.01753715, 0.01833317,0.01814674, 0.01594022, 0.01273109, 0.01097039, 0.01016107,0.00947185, 0.00875745]]))

# Geometries of view and illumination--------------------------------------------------------------------------:
s.geometry = Geometry.User()
s.geometry.day = 3
s.geometry.month = 10
s.geometry.solar_z = float(solar_Zn)
s.geometry.solar_a = float(solar_Az)
s.geometry.view_z = float(view_Zn)

# Altitudes---------------------------------------------------------------------------------------------------:
s.altitudes = Altitudes()
s.altitudes.set_sensor_satellite_level()  # Set the sensor altitude to be satellite level.
s.altitudes.set_target_custom_altitude(alt_alvo)  # The altitude of the target (in km).

#%%
# Função de interpolação 
def _interpolate_(dataframe_):
    import scipy.interpolate
    a = np.arange(379,1049,2.5) #Para o VNIR
    #a = np.arange(907,2519,2.5) #Para o SWIR
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
# Função de resposta espectral do PRISMA
frs = pd.read_excel('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA Spectral.xlsx',index_col=0)

frs_nova = _interpolate_(frs)
frs_nova_transposed = frs_nova.transpose()
frs = frs.transpose()

#%%
# Rodar o 6S para cada banda espectral    
dicionario = {}
waves_6s = []
waves_max_int = []

for a in range(0,frs_nova_transposed.shape[1]):
    
    max_value_index = frs_nova_transposed[frs_nova_transposed.iloc[:,a]==frs_nova_transposed.iloc[:,a].max()].index.values[0]
    waves_max_int.append(max_value_index)
    
    wv_value = frs_nova_transposed.columns.values[a]
    
    # Comparação dos comprimentos de onda (não usa pra nada)
    max_value_temp = frs[frs.iloc[:,a]==frs.iloc[:,a].max()].index.values[0]
    waves_6s.append(max_value_temp)
    
    #wv_init = max_value_index - 22.5
    
    #vw_final = max_value_index + 22.5
    
    srf_values = list(frs_nova_transposed.loc[379:1046.5,wv_value])
    #srf_values = list(frs_nova_transposed.loc[907:2517.5,wv_value])
    
    s.wavelength = Wavelength((379/1000),(1046.5/1000),srf_values) # Para o VNIR
    #s.wavelength = Wavelength((907/1000),(2517.5/1000),srf_values) # Para o SWIR
    
    s.run()
    
    print(s.wavelength)
    
    # Alguns parâmetros importantes não estão contidos no output do modelo. Podemos recuperar esses valores:
    #s.outputs.values['transmittance_total_scattering'] = s.outputs.transmittance_total_scattering.total
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
    
    dicionario[str(int(max_value_temp))] = s.outputs.values #Guarda dentro de um dicionário, para cada banda, os dados de saída do 6S

#%%
# Exportar txt com todas as variáveis possíveis
#s.outputs.write_output_file('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/py6s_output.txt')

#%%
# Salvar dicionário em csv
#with open('G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/PRISMA_04-09-21_resultados/6S_parameters_SWIR.csv', 'w') as f:
#    for key in dicionario.keys():
#        f.write("%s,%s\n"%(key,dicionario[key]))
    
#%%
# Salva dicionário em json
import json
a_file = open("G:/Outros computadores/Meu modelo Laptop (1)/Documents/Mestrado/BEPE/atmospheric_correction/6S/testes_setembro/6S_parameters_0.060.json", "w")
json.dump(dicionario, a_file)
