# Python 3.5
# Script for MetalSilicate python library
# Copyright Kayla Iacovino
# This script contains pandas dataframes with interaction parameter data digitized from the Steelmaking Data Source Book

import pandas as pd

# These hard-coded pandas DataFrames are generated from an Excel Spreadsheet and hardcoded for easier access.
"""
The code to generate these:
import pandas as pd
import numpy as np
interaction_parameter_df = pd.read_excel('ActivityCoeffs_and_InteractionParams.xlsx', sheet_name='InteractionParameters_Liquid')
interaction_parameter_df.fillna("No Data")

import math
f= open("hard_coded_calibrations.py","w+")
cols = ['i,j,k', 'ei(j,k)', 'Temp_K', 'ConcRange_MassPercent', 'TempDependency', 'TempRange_K', 'Reference', 'Note']

f.write("\n")
f.write("interaction_params = pd.DataFrame({ \n")
for col in cols:
	iterno = 1
	f.write("'" + (str(col)+"': ["))
	for index, row in interaction_parameter_df.iterrows():
		try:
			value = float(row[col])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[col]) + "'")
		if iterno < len(interaction_parameter_df.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }).set_index('i,j,k') \n")


##NOTE! Delete the final comma at the end of the final column.

##Reminder: to look up a value, use syntax: interaction_params.loc['Ag Ag']['Temp_K']
"""


interaction_params = pd.DataFrame({ 
'i,j,k': ['Ag Ag','Ag Al','Ag C','Ag Cr','Ag O','Al Ag','Al Al','Al C','Al Ca','Al H','Al N','Al O','Al P','Al Pb','Al S','Al Si','Al U','As C','As N','As S','Au O','Au S','B B','B C','B H','B Mn','B N','B O','B P','B S','B Si','Be O','C Ag','C Al','C As','C B','C C','C Ca','C Ce','C Co','C Co,C','C Co,Co','C Cr','C Cu','C Ge','C H','C La','C Mg','C Mn','C Mo','C Mo,Mo','C N','C Nb','C Nb,Nb','C Ni','C Ni,C','C Ni,Ni','C O','C P','C Pb','C S','C Sb','C Si','C Sn','C Ta','C V','C V,C','C V,V','C W','C W,W','Ca Al','Ca C','Ca Ca','Ca Ni','Ca O','Ca S','Ca Si','Ce C','Ce Ce','Ce H','Ce Mn','Ce O','Ce S','Co C','Co Co','Co Co,Co','Co Cr','Co H','Co Mn','Co N','Co O','Co P','Co Pb','Co S','Cr Ag','Cr C','Cr Co','Cr Cr','Cr Cu','Cr H','Cr Mn','Cr Mo','Cr N','Cr Ni','Cr O','Cr P','Cr Pb','Cr S','Cr Si','Cr Sn','Cr V','Cu C','Cu Cr','Cu Cu','Cu Cu,Cu','Cu H','Cu N','Cu O','Cu P','Cu Pb','Cu S','Cu Si','Cu Sn','Ge C','Ge Ge','Ge H','Ge S','H Al','H B','H C','H Ce','H Co','H Cr','H Cu','H Ge','H H','H La','H Mn','H Mo','H Nb','H Nd','H Ni','H O','H P','H Pd','H Rh','H S','H Si','H Sn','H Ta','H Ti','H V','H W','H Zr','Hf Hf','Hf O','Hf S','La C','La H','La La','La Mn','La O','La S','Mg C','Mg O','Mg B','Mn C','Mn Ce','Mn Co','Mn Cr','Mn H','Mn La','Mn Mn','Mn Mo','Mn N','Mbn Nb','Mn Ni','Mn O','Mn P','Mn Pb','Mn S','Mn Si','Mn Ta','Mn Ti','Mn V','Mn W','Mo C','Mo Cr','Mo H','Mo Mn','Mo Mo','Mo Mo,Mo','Mo N','Mo O','Mo P','Mo Pb','Mo S','Mo Si','N Al','N As','N B','N C','N Co','N Cr','N Cr,Cr','N Cu','N Mn','N Mo','N Mo,Mo','N N','N Nb','N Ni','N O','N P','N S','N Sb','N Se','N Si','N Sn','N Ta','N Ta,Ta','N Ta','N Ta,Ta','N Te','N Ti','N V','N V','N V,V','N W','N Zr','Nb C','Nb H','Nb Mn','Nb N','Nb Nb','Nb O','Nb P','Nb S','Nb Si','Ni C','Ni Ca','Ni Cr','Ni H','Ni Mn','Ni N','Ni Ni','Ni O','Ni P','Ni Pb','Ni S','Ni Si','O Ag','O Al','O Au','o B','O Be','O C','O Ca','O Ce','O Co','O Cr','O Cu','O H','O Hf','O La','O Mg','O Mn','O Mo','O N','O Nb','O Ni','O O','O P','O Pd','O Pt','O Rh','O S','O Sb','O Sc','O Si','O Sn','O Ta','O Ti','O U','O V','O W','O Y','O Zr','P Al','P B','P C','P Co','P Cr','P Cu','P H','P Mn','P Mo','P N','P Nb','P Ni','P O','P P','P Pb','P S','P Si','P Sn','P Ti','P V','P W','Pb Al','Pb C','Pb Co','Pb Cr','Pb Cu','Pb Mn','Pb Mo','Pb NI','Pb P','Pb S','Pb Si','Pb Sn','Pb W','Pd H','Pd O','Pd Pd','Pt O','Pt S','Rh H','Rh O','S Al','S As','S Au','S B','S C','S Ca','S Ce','S Co','S Cr','S Cu','S Ge','S H','S Hf','S La','S Mn','S Mo','S N','S Nb','S Ni','S O','S P','S Pb','S Pt','S S','S Sb','S Si','S Sn','S Ta','S Ti','S U','S V','S W','S Y','S Zr','Sb C','Sb N','Sb O','Sb S','Sc O','Se N','Si Al','Si B','Si C','Si Ca','Si Cr','Si Cu','Si H','Si Mn','Si Mo','Si N','Si Nb','Si Ni','Si O','Si P','Si Pb','Si S','Si Si','Si Sn','Si Ta','Si Ti','Si V','Sn C','Sn Cr','Sn Cu','Sn H','Sn N','Sn O','Sn P','Sn Pb','Sn S','Sn Si','Sn Sn','Ta C','Ta H','Ta Mn','Ta N','Ta O','Ta S','Ta Si','Ta Ta','Te N','Ti H','Ti Mn','Ti N','Ti O','Ti P','Ti S','Ti Si','Ti Ti','U Al','U O','U S','U U','V C','V Cr','V H','V Mn','V N','V N','V O','V P','V S','V Si','V V','W C','W H','W Mn','W N','W O','W P','W Pb','W S','Y O','Y S','Zr H','Zr N','Zr O','Zr S','Zr Zr'], 
'ei(j,k)': [-0.4,-0.08,0.22,-0.0097,-0.099,-0.017,0.043,0.091,-0.047,0.24,0.015,-1.98,0.033,0.0065,0.035,0.056,0.011,0.25,0.077,0.0037,-0.14,-0.0051,0.038,0.22,0.58,-0.00086,0.073,-0.21,0.008,0.048,0.078,-1.3,0.028,0.043,0.043,0.244,0.243,-0.097,-0.0026,0.0075,0.0002,0.0009,-0.023,0.016,0.008,0.67,0.0066,0.07,-0.0084,-0.0137,-6e-05,0.11,-0.059,-0.00024,0.01,0.00026,0.0015,-0.32,0.051,0.0099,0.044,0.0154,0.08,0.022,-0.23,-0.03,0.00029,0.012,-0.0056,-5e-05,-0.072,-0.34,-0.002,-0.044,-580.0,-140.0,-0.096,-0.077,0.0039,-0.6,0.13,-560.0,-40.0,0.02,0.00509,2.4e-05,-0.022,-0.14,-0.0042,0.037,0.018,0.0037,0.0031,0.0011,-0.0024,-0.114,-0.019,-0.0003,0.016,-0.34,0.0039,0.0018,-0.182,0.0002,-0.16,-0.033,0.0083,-0.17,-0.004,0.009,0.012,0.066,0.018,-0.02,0.00013,-0.19,0.025,-0.065,-0.076,-0.0056,-0.021,0.027,-0.011,0.03,0.007,0.41,0.026,0.013,0.058,0.06,0.0,0.0018,-0.0024,0.0013,0.01,0.0,-0.027,-0.002,0.0029,-0.0033,-0.038,-0.0019,0.05,0.015,0.0041,0.056,0.017,0.027,0.0057,0.0017,-0.019,-0.0074,0.0048,-0.0088,0.007,-3.2,-0.27,0.03,-4.3,-0.0078,0.28,-43.0,-79.0,0.15,-3.0,-0.0236,-0.0538,0.054,-0.0036,0.0039,-0.34,0.11,0.0,0.0046,-0.091,0.0073,-0.0072,-0.083,-0.06,-0.0029,-0.048,-0.0327,0.0035,-0.05,0.0057,0.0071,-0.14,-0.0003,-0.13,0.0048,0.0121,-0.00011,-0.1,0.0083,-0.006,0.0023,-0.0006,8.05,0.01,0.018,0.094,0.13,0.012,-0.046,0.00028,0.009,-0.02,-0.011,-0.00048,0.0,-0.068,0.007,-0.12,0.059,0.007,0.0088,0.006,0.048,0.007,-0.049,0.0,-0.058,0.0082,0.07,-0.6,-0.123,-0.111,0.002,-0.002,-0.63,-0.486,-0.7,0.0093,-0.475,0.0,-0.72,-0.045,-0.046,-0.01,0.032,-0.066,-0.0003,-0.36,-0.008,0.015,0.0007,0.01,0.0018,-0.0023,-0.0036,0.006,-0.011,-1.17,-0.007,-0.31,-2.4,-0.421,-515.0,-64.0,0.008,-0.055,-0.013,0.73,-0.28,'0-5',-1.98,-0.021,0.005,-0.14,-0.12,0.006,-0.17,0.07,-0.009,0.0045,0.0136,-0.133,-0.023,-1.3,-0.066,-0.0111,-0.1,-1.12,-0.44,-0.14,0.0085,-0.46,-4.0,0.037,0.015,0.126,0.04,-0.018,-0.035,0.33,-0.032,0.001,0.13,-0.012,0.003,0.13,0.054,0.011,0.034,0.099,0.013,-0.04,-0.024,-0.023,0.021,0.1,0.0,0.02,-0.028,-0.023,0.0,-0.019,0.048,-0.32,0.048,0.057,0.0,0.021,-0.084,0.002,0.0063,0.032,0.13,0.064,0.041,0.0041,0.0028,0.134,0.111,-110.0,-9.1,0.0026,-0.0105,-0.0084,0.014,0.41,-0.045,-18.3,-0.026,0.0027,0.01,-6.013,0.0,-0.27,0.035,-0.046,0.0089,-0.046,0.0037,0.075,-0.0044,-0.019,-0.18,-0.067,-0.019,0.011,-0.275,-0.21,0.11,0.043,-0.2,0.0019,-3.7,0.014,0.058,0.2,0.18,-0.066,-0.0003,0.0144,0.64,-0.0146,2.36,0.092,0.0,0.005,-0.119,0.09,-0.01,0.066,0.103,0.017,0.04,1.23,0.025,0.18,0.015,-0.024,0.16,0.027,-0.11,0.036,0.035,-0.028,0.057,0.0017,-3.5,-0.47,0.0016,-0.685,-1.2,-0.13,0.23,0.11,0.6,-1.1,-0.043,-2.06,-3.4,-0.06,-0.27,-0.27,0.042,0.059,-6.6,-0.53,0.013,-0.14,0.0119,-0.59,0.0056,-0.455,-0.4,-0.46,-0.042,-0.033,0.042,0.0309,-0.15,0.088,0.0136,-0.079,0.052,-0.16,0.0005,0.043,-2.6,-0.77,-1.2,-4.13,-23.0,-0.61,0.032], 
'Temp_K': [1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1880.0,1873.0,1953.0,1873.0,1873.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1823.0,1823.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1673.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1823.0,1823.0,1873.0,1833.0,1873.0,1873.0,1873.0,1873.0,1843.0,1833.0,1833.0,1873.0,1833.0,1833.0,1823.0,1823.0,1823.0,1873.0,1873.0,1823.0,1823.0,1873.0,1873.0,1823.0,1833.0,1823.0,1823.0,1823.0,1833.0,1833.0,1880.0,1873.0,1873.0,1880.0,1873.0,1873.0,1880.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1873.0,1903.0,1873.0,1843.0,1873.0,1873.0,1873.0,1823.0,1823.0,1873.0,1873.0,1903.0,1873.0,1873.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1903.0,1823.0,1873.0,1833.0,1873.0,1873.0,1873.0,1873.0,1879.0,1873.0,1673.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1953.0,1883.0,1873.0,1873.0,1873.0,1873.0,1873.0,1843.0,1843.0,1873.0,1873.0,1863.0,1843.0,1873.0,1843.0,1843.0,1873.0,1673.0,1823.0,1823.0,1873.0,1834.0,1873.0,1843.0,1843.0,1833.0,1873.0,1873.0,1843.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1863.0,1873.0,1953.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1879.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1853.0,1873.0,1853.0,1873.0,1879.0,1873.0,1873.0,1873.0,1873.0,1853.0,1873.0,1873.0,1873.0,1873.0,1879.0,1873.0,1833.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1823.0,1880.0,1873.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1823.0,1823.0,1873.0,1873.0,1873.0,1823.0,1873.0,1873.0,'<1.0',1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1953.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1874.0,1873.0,1873.0,1873.0,1873.0,1873.0,1673.0,1873.0,1873.0,1873.0,1673.0,1873.0,1673.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1723.0,1823.0,1873.0,1823.0,1873.0,1873.0,1673.0,1823.0,1823.0,1823.0,1873.0,1823.0,1823.0,1873.0,1823.0,1723.0,1723.0,1823.0,1823.0,1823.0,1873.0,1823.0,1873.0,1873.0,1823.0,1873.0,1823.0,1823.0,1823.0,1823.0,1823.0,1823.0,1873.0,1873.0,1823.0,1823.0,1823.0,18213.0,1873.0,1873.0,1883.0,1823.0,1873.0,1853.0,1823.0,1823.0,1873.0,1823.0,1723.0,1823.0,1873.0,1823.0,1823.0,1823.0,1873.0,1873.0,1873.0,1823.0,1825.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1853.0,1873.0,1873.0,1873.0,1880.0,1903.0,1873.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1823.0,1873.0,1823.0,1873.0,1873.0,1833.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1823.0,1823.0,1823.0,1823.0,1873.0,1833.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1873.0,1853.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1823.0,1873.0,1873.0,1843.0,1873.0,1873.0,1873.0,1873.0,1823.0,1833.0,1873.0,1833.0,1873.0,1843.0,1879.0,1873.0,1673.0,1823.0,1823.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0], 
'ConcRange_MassPercent': ['No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','Al < 3.8','No Data','No Data','No Data','Al < 7','No Data','No Data','No Data','No Data','No Data','No Data','No Data','< 7.06','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<20','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<11','No Data','No Data','<20','No Data','No Data','<10.5','<3.4','<2.0','No Data','<2.1','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','Cr<60','No Data','No Data','No Data','No Data','No Data','<5','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<4','<2.5','No Data','<0.6','<14','<15','No Data','<11','No Data','<1.6','No Data','No Data','<15','No Data','No Data','No Data','<6','No Data','No Data','<1','<10','<7','<25','<2','<20','<20','<2','No Data','No Data','No Data','<0.7','La<1.6','<0.1','<1.0','No Data','No Data','No Data','No Data','<2.8','<1.2','No Data','<4.3','<3.2','No Data','No Data','No Data','<3.2','No Data','<4.5','<3.5','No Data','No Data','No Data','No Data','<2.8','<4.6','No Data','<3.5','<3.3','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<3.8','<2','<7.06','<2','<5.5','<60','<60','<9','<4','<8','No Data','No Data','<9.7','<5.5',0.12,'<1.1','<4','<5','<4','<3.5','<4','<7.1','<7.1','7.1~20','7.1~20','<0.6','<0.5','<4.8','4.8~17.9','4.8~17.9','<15','<0.6','No Data','No Data','No Data','Nb<9.7','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','0.32~2.97','No Data','No Data','No Data','No Data','<40','No Data','<5','No Data','0.07~1.3','<0.6','No Data','No Data','<10','No Data','<3','<40','No Data','<0.7','<15','<30','<20','No Data','<14','No Data','<3','<14','<2','<1','0.003~4.87','<12','<22','0.01~1.85','<0.15','No Data','<3.7','No Data','No Data','No Data','<9.9','P<6','<19.3','No Data','No Data','No Data','No Data','<0.7','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<15','No Data','<20','No Data','No Data','<5.5','No Data','No Data','No Data','No Data','<15','<15','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<7','No Data','No Data','No Data','No Data','No Data','<3','No Data','No Data','No Data','No Data','No Data','1.8~16','<0.1','No Data','zMo<0.3','No Data','zNb<0.1','No Data','S<1.2','<9','No Data','No Data','No Data','No Data','<7','No Data','<2.5','1.8~6','0.17~17','<11','<15','0.2-8.35','0.2-3.0','No Data','No Data','Sb~14','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','Si<3','No Data','No Data','Si<7','<3','No Data','No Data','No Data','<4.5','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','<4.6','<7.1','No Data','No Data','No Data','<2','No Data','No Data','<2.2','No Data','No Data','No Data','No Data','No Data',0.023,'0.4~4','No Data','No Data','<0.8','No Data','No Data','No Data','No Data','<4.8','4.8~17.9','No Data','No Data','V<11','No Data','<24','No Data','No Data','No Data','No Data','No Data','No Data','No Data','W<15','No Data','No Data','No Data','No Data','No Data','No Data','No Data'], 
'TempDependency': ['No Data','No Data','No Data','No Data','No Data','No Data','80.5/T','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','700/T - 0.337','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','162/T-0.008','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-549/T+0.111','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-86.9/T+0.0336','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-37.3/T+0.0166','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-126/T+0.0458','No Data','No Data','-76.5/T+0.0321','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','-1370/T+0.690','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','413/T-0.217','No Data','No Data','No Data','No Data','No Data','-1838/T+0.964','No Data','No Data','No Data','236/T-0.120','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','-148/T+0.033','1.56/T-0.00053','No Data','No Data','-33.2/T+0.0064','-5.57/T+0.0025','No Data','-280/T+0.0816','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-524/T+0.231','No Data','-702/T+0.317','21.2/T-0.0104','No Data','-5700/T+2.45','-1420/T+0.635','-1220/T+0.544','30/T-0.014','No Data','*','No Data','No Data','No Data','-1860/T+0.517','No Data','-19970/T+9.950','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-3440/T+1.717','No Data','-1750/T+0.76','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-1830/T+0.874','No Data','No Data','-1050/T+0.42','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','-94.2/T+0.040','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-120/T+0.018','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','380/T-0.023','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','-162.4/T+0.0627','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-6770/T+2.93','-20700/T+9.83','No Data','No Data','4737/T-2.42','No Data','No Data','No Data','-19500/T+8.37','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','-5160/T+2.3','-4440/T+1.97','No Data','No Data','No Data','No Data','30/T-0.22','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','*','No Data'], 
'TempRange_K': ['No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1733–2033','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','*','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1833–1883','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1843–1868','No Data','No Data','No Data','No Data','No Data','1843–1868','No Data','No Data','1821–1945','1843–1868','No Data','No Data','1700–1883','No Data','No Data','1823–1873','1813–1878','1821–1945','1821–1945','1821–1945','No Data','1821–1945','1821–1945','No Data','No Data','No Data','No Data','1843–1868','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1823–1873','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1823–1873','No Data','No Data','No Data','No Data','No Data','1818–1893','No Data','No Data','No Data','1823–1873','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1823–2053','1823–2053','No Data','1823–1973','1773–2073','1773–2073','No Data','1873–1973','1873–1973','1853–1953','No Data','No Data','No Data','No Data','1823–1973','No Data','1873–1973','No Data','1873–1973','1873–1973','No Data','1873–1973','1873–1973','1873–1973','1873–1973','No Data','No Data','No Data','No Data','No Data','1873–1973','1873–1973','1823–1923','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1733–2033','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1823–1923','No Data','1823–1923','No Data','No Data','1843–1873','No Data','No Data','No Data','No Data','1823–1973','No Data','No Data','No Data','No Data','1823–1923','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1700–1883','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1873–1923','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1818–1893','No Data','No Data','No Data','No Data','1823–1973','No Data','No Data','No Data','1823–1973','No Data','No Data','No Data','No Data','No Data','No Data','1833–1883','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1873–1973','No Data','No Data','No Data','1823–1923','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','1873–1973','1873–1973','1823–1923','No Data','No Data','No Data','No Data','No Data','No Data','1823–1873','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data'], 
'Reference': [1.0,4.0,2.0,1.0,5.0,4.0,13.0,10.0,17.0,20.0,31.0,13.0,39.0,41.0,46.0,10.0,2.0,49.0,50.0,43.0,5.0,43.0,53.0,49.0,21.0,57.0,53.0,35.0,61.0,43.0,2.0,62.0,2.0,10.0,49.0,49.0,90.0,2.0,102.0,105.0,105.0,105.0,121.0,86.0,2.0,130.0,102.0,49.0,134.0,94.0,94.0,138.0,94.0,94.0,105.0,105.0,105.0,90.0,1.0,40.0,44.0,171.0,2.0,92.0,94.0,105.0,105.0,105.0,94.0,94.0,17.0,2.0,2.0,17.0,193.0,192.0,17.0,102.0,102.0,195.0,102.0,199.0,102.0,105.0,203.0,203.0,204.0,22.0,134.0,212.0,214.0,39.0,41.0,43.0,1.0,121.0,204.0,2.0,124.0,205.0,57.0,2.0,237.0,2.0,252.0,39.0,2.0,259.0,223.0,124.0,260.0,86.0,124.0,263.0,263.0,208.0,25.0,216.0,266.0,41.0,43.0,124.0,513.0,2.0,268.0,22.0,43.0,22.0,21.0,130.0,195.0,22.0,205.0,208.0,22.0,270.0,195.0,206.0,205.0,207.0,195.0,20.0,206.0,129.0,208.0,208.0,19.0,18.0,21.0,21.0,207.0,205.0,207.0,21.0,281.0,201.0,201.0,102.0,195.0,102.0,102.0,198.0,102.0,49.0,62.0,57.0,154.0,102.0,134.0,57.0,206.0,102.0,284.0,57.0,141.0,57.0,134.0,2.0,266.0,41.0,43.0,134.0,57.0,102.0,134.0,57.0,94.0,2.0,205.0,57.0,298.0,298.0,293.0,214.0,39.0,2.0,43.0,47.0,31.0,50.0,53.0,138.0,212.0,237.0,237.0,25.0,141.0,293.0,209.0,306.0,238.0,210.0,143.0,211.0,322.0,50.0,323.0,141.0,25.0,238.0,238.0,238.0,238.0,323.0,329.0,238.0,238.0,238.0,25.0,230.0,94.0,207.0,57.0,238.0,311.0,342.0,39.0,43.0,427.0,105.0,17.0,2.0,20.0,134.0,210.0,202.0,215.0,39.0,41.0,46.0,151.0,5.0,13.0,5.0,35.0,62.0,90.0,199.0,199.0,214.0,252.0,216.0,206.0,201.0,198.0,62.0,2.0,214.0,143.0,342.0,215.0,360.0,368.0,370.0,51.0,370.0,91.0,243.0,374.0,380.0,255.0,342.0,7.0,374.0,392.0,51.0,201.0,188.0,39.0,61.0,1.0,39.0,39.0,266.0,129.0,266.0,39.0,211.0,39.0,39.0,368.0,401.0,41.0,46.0,39.0,124.0,39.0,39.0,266.0,41.0,40.0,41.0,2.0,41.0,41.0,2.0,41.0,41.0,41.0,41.0,41.0,41.0,208.0,370.0,403.0,51.0,43.0,208.0,370.0,46.0,43.0,43.0,43.0,44.0,192.0,102.0,43.0,259.0,43.0,43.0,19.0,201.0,102.0,43.0,43.0,322.0,43.0,46.0,91.0,46.0,41.0,43.0,412.0,43.0,46.0,124.0,201.0,102.0,201.0,46.0,46.0,201.0,102.0,171.0,50.0,243.0,43.0,374.0,323.0,10.0,2.0,2.0,17.0,223.0,124.0,18.0,134.0,47.0,141.0,427.0,151.0,380.0,39.0,41.0,46.0,380.0,124.0,427.0,47.0,423.0,92.0,124.0,513.0,21.0,25.0,255.0,124.0,41.0,124.0,124.0,513.0,94.0,21.0,57.0,238.0,342.0,201.0,427.0,342.0,323.0,207.0,102.0,329.0,7.0,39.0,102.0,47.0,390.0,48.0,374.0,201.0,2.0,105.0,260.0,205.0,134.0,238.0,238.0,392.0,39.0,46.0,423.0,428.0,94.0,207.0,57.0,25.0,51.0,266.0,41.0,46.0,201.0,201.0,21.0,230.0,188.0,102.0,32.0], 
'Note': ['assum., quasi-chem. [C]','from eAl(Ag) [B]','[C]','from eCr(Ag) [C]','from eO(Ag) [C]','[B]','[A]','from eAl(C) [C]','from eCa(Al) [B]','from eH(Al) [A]','from eN(Al) [B]','from eO(Al) [A]','from eP(Al) [A]','[B]','form eS(Al) [A]','from e [A]','from eU(Al) [C]','from eC(As) [C]','from eN(As) [B]','from eS(As) [A]','from eO(Au) [A]','from eS(Au) [A]','calcd. [C]','from eC(B) [C]','from eH(B) [A]','from eMn(B) [C]','from eN(B) [C]','from eO(B) [B]','from thetaP(B) [B]','from eS(B) [A]','[C]','from eO(Be) [C]','[C]','from eAl(C) [C]','[C]','calcd. [C]','gamma_Prime_C given by func. [A]','[C]','from eCe(C) [C]','from e [B]','conv. [B]','conv. [B]','[A]','[B]','[C]','from eH(C) [A]','from eLa(C) [C]','from eC(Mg) [C]','from Mn(C) [B]','from e [B]','conv. [B]','from eN(C) [A]','from e [B]','conv [B]','from e [A]','conv. [A]','conv. [A]','from eO(C) [A]','from eP(C) [A]','from ePb(C) [B]','from eS(C) [A]','[C]','[A]','from e [B]','from e [C]','[B]','[B]','[B]','[B]','[B]','from eCa(Al) [B]','[C]','[B]','from eCa(Ni) [B]','from eO(Ca) [B]','from eS(Ca) [A]','from eCa(Si) [A]','[C]','from eCe(Ce) [A]','from eH(Ce) [C]','[B]','from eO(Ce) [C]','from eS(Ce) [B]','from eC(Co) [B]','from e gammaCo given [A]','from given gammaCo [A]','from eCo(Cr) [C]','from eH(Co) [A]','from eMn(Co) [A]','from eN(Co) [A]','from eO(Co) [A]','from eP(Co) [A]','from ePb(Co) [B]','from eS(Co) [A]','from e [C]','from eC(Cr) [A]','from eCr(Co) [C]','[A]','from eCr(Cu) [C]','from eH(Cr) [A]','from eMn(Cr) [A]','[C]','from eN(Cr) [A]','[C]','from eO(Cr) [A]','from eP(Cr) [A]','[C]','from eS(Cr) [A]','[B]','from eCr(Sn) [C]','[B]','from eC(Cu) [B]','from eCu(Cr) [C]','[A]','[A]','from eH(Cu) [A]','from eN(Cu) [A]','from eO(Cu) [A]','from eP(Cu) [B]','from ePb(Cu) [B]','from eS(Cu) [A]','from eCu(Si) [B]','from eSn(Cu) [B]','[C]','[C]','from eH(Ge) [A]','from thetaS(Ge) [A]','from Fig 2 [A]','[A]','[A]','calcd. In ref 2, [C]','from Fig. 2 [A]','[A]','from e [A]','from Fig. 2 [A]','[A]','from Fig. [C]','calcd. By eq. [A]','[A]','[A]','from Fig. [C]','from e [A]','from solubility [C]','[A]','from e [A]','from e [A]','[A]','[A]','[A]','[B]','[A]','[A]','[A]','[A]','from eHf(Hf) [C]','from eO(Hf) [C]','from eS(Hf) [C]','[C]','from eH(La) [C]','from eLa(La) [B]','[B]','from eO(La) [C]','from eS(La) [A]','from eC(Mg) [C]','from eO(Mg) [C]','[C]','[B]','from eCe(Mn) [B]','[A]','[A]','from eH(Mn) [A]','from eLa(Mn) [B]','conv. [A]','[A]','from eN(Mn) [A]','±0.003 [B]','[A]','[A]','from eP(Mn) [B]','from ePb(Mn)','from eS(Mn) [A]','[A]','±0.003 [B]','from eTi(Mn) [C]','[A]','[A]','from eC(Mo) [B]','[C]','from eH(Mo) [A]','from eMn(Mo) [A]','±0.0035 [C]','±0.00008 [C]','from eN(Mo) [A]','from eO(Mo) [A]','from eP(Mo) [A]','[C]','from eS(Mo) [A]','[C]','[B]','[B]','[C]','[A]','bN(Co) [A]','from e [A]','conv. [A]','[A]','[A]','[A]','[A]','[A]','[C]','[A]','[B]','[A]','[B]','[B]','[B]','[A]','[B]','[C]','[C]','[C]','[C]','[B]','[B]','[A]','[A]','[A]','[A]','[C]','from eC(Nb) [B]','from eH(Nb) [A]','from eMn(Nb) [B]','from eN(Nb) [C]','calcd. Erroe ±0.009 [A]','from eO(Nb) [A]','from eP(Nb) [A]','from eS(Nb) [A]','from eSi(Nb) [B]','from eC(Ni) [A]','from eCa(Ni) [B]','[C]','from eH(Ni) [A]','from eMn(Ni) [A]','from eN(Ni) [A]','from eNi(Ni) [A]','from eO(Ni) [A]','from eP(Ni) [A]','from ePb(Ni) [B]','from eS(Ni) [A]','from eSi(Ni) [B]','[C]','[A]','±0.0005 [A]','[B]','[C]','[A]','[A]','[B]','[A]','[A]','[A]','from eH(O) [C]','[C]','EMF [C]','[C]','[A]','[A]','from eN(O) [B]','[A]','[A]','[A]','[B]','from Fig. 3 [B]','from Fig. 6 [A]','from Fig. 2 [B]','from eS(O) [A]','[B]','[C]','[A]','[B]','[A]','EMF [A]','[C]','[A]','[A]','[C]','[B]','from eP(Al) [A]','from thetaP(B) [B]','[A]','from eP(Co) [A]','from eP(Cr) [A]','from thetaP(Cu) [B]','from eH(P) [A]','from thetaP(Mn) [B]','from eP(Mo) [A]','from eN(P) [A]','from eP(Nb) [A]','from eP(Ni) [A]','from eO(P) [B]','from eP(P) [A]','from ePb(P) [B]','from eS(P) [A]','from eP(Si) [A]','from eSN(P) [B]','from eP(Ti) [A]','from eP(V) [A]','[B]','[B]','[B]','[B]','[C]','±0.002 [B]','[B]','[C]','±0.001 [B]','[B]','[B]','[B]','[B]','[B]','from eH(Pd) [A]','from eO(Pd) [B]','estimated [C]','from eO(Pt) [A]','from eS(Pt) [A]','from eH(Rh) [A]','from eO(Rh) [B]','[A]','from thetaS(As) [A]','from thetaS(Au) [A]','from thetaS(B) [A]','[A]','calcd. Value [A]','[A]','from thetaS(Co) [A]','[A]','from thetaS(Cu) [A]','from thetaS(Ge) [A]','from eH(S) [A]','[C]','[A]','from thetaS(Mn) [A]','from thetaS(Mo) [A]','from eN(S) [B]','[A]','[A]','[A]','[A]','from ePb(S)','from thetaS(Pt) [A]','[A]','from thetaS(Sb) [A]','[A]','from eSn(S) [A]','[B]','[B]','[C]','[A]','[A]','[C]','[B]','from eC(Sb) [C]','from eN(Sb) [B]','from eO(Sb) [B]','from eS(Sb) [A]','from eO(Sc) [C]','from eN(Se) [B]','from eAl(Si) [A]','[C]','[A]','from eCa(Si) [A]','from eCr(Si) [B]','from eCu(Si) [B]','from eH(Si) [A]','from eMn(Si) [A]','[C]','from eN(Si) [A]','[B]','[B]','[A]','from eP(Si) [A]','from ePb(Si) [B]','from eS(Si) [A]','[A]','from eSn(Si) [B]','[B]','[C]','from e [A]','from eC(Sn) [B]','from eSn(Cr) [C]','from eSn(Cu) [B]','from eH(Sn) [A]','from eN(Sn) [B]','from eO(Sn) [B]','from eP(Sn) [B]','from ePb(Sn) [B]','from eSn(s) [A]','from e [B]','from e Mass-spec [B]','from eC(Ta) [C]','from eH(Ta) [B]','from eMn(Ta) [B]','from eN(Ta) [C]','from eO(Ta) [A]','from eS(Ta) [B]','from eSi(Ta) [B]','from eq. [A]','from eN(Te) [B]','from eH(Ti) [A]','[C]','from eN(Ti) [B]','from eO(Ti) [A]','from eP(Ti) [A]','from eS(Ti) [B]','[C]','from e [A]','from e [C]','from eO(U) [C]','from eS(U) [B]','from e [B]','from eC(V) [B]','extrapolated to V=0, Cr=0 [B]','from eH(V) [A]','from eMn(V) [A]','from eN(V) [A]','from eN(V) [A]','from eO(V) [A]','from eP(V) [A]','from eS(V) [A]','from eSi(V) [A]','[A]','from eC(W) [B]','from eH(W) [A]','from eMn(W) [A]','from eN(W) [A]','from eO(W) [A]','from eP(W) [B]','from ePb(W) [B]','from eS(W) [A]','from eO(Y) [B]','from eS(Y) [B]','from eH(Zr) [A]','from eN(Zr) [C]','from eO(Zr) [B]','from eS(Zr) [B]','EMF [C]'] 
 }).set_index('i,j,k')

