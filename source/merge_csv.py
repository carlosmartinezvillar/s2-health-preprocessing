'''
File to merge and filter CSV files from CDC data and USDA RUCA.
'''

'''
Keep columns:
from RUCA: 
	TractFIPS20, 
	CountyName20, 
	StateName20,
	PrimaryRUCA, 
	Population, 
	LandArea (float), 
	PopDensity (float)

from CDC: 
	StateDesc (str), 
	CountyName (str), 
	LocationID (as str),
	Data_Value (float), 
	Low_Confidence_Limit (float), 
	High_Confidence_Limit (float),
	TotalPopulation (int), 
	TotalPop18plus (int)

Merged by TractFIPS20-LocationID	
'''
# LIBRARIES
import pandas as pd

# KEEP ONLY THESE COLUMNS
ruca_cols = ['TractFIPS20', 'CountyName20', 'StateName20',
             'PrimaryRUCA', 'Population', 'LandArea', 'PopDensity']
cdc_cols  = ['StateDesc', 'CountyName', 'LocationID',
             'Data_Value', 'Low_Confidence_Limit', 'High_Confidence_Limit',
             'TotalPopulation', 'TotalPop18plus']

# LOAD FILES
with open('../shapes/RUCA-codes-2020-tract.csv',encoding='utf-8',errors='replace') as fp:
	ruca = pd.read_csv(fp,usecols=ruca_cols)
cdc  = pd.read_csv('../shapes/cdc_diabetes.csv', usecols=cdc_cols)

# CAST TYPES
ruca['LandArea']   = ruca['LandArea'].astype(float)
ruca['PopDensity'] = ruca['PopDensity'].astype(float)
cdc['Data_Value']             = cdc['Data_Value'].astype(float)
cdc['Low_Confidence_Limit']   = cdc['Low_Confidence_Limit'].astype(float)
cdc['High_Confidence_Limit']  = cdc['High_Confidence_Limit'].astype(float)
cdc['TotalPopulation']        = cdc['TotalPopulation'].apply(lambda s: s.replace(",",""))
cdc['TotalPop18plus']         = cdc['TotalPop18plus'].apply(lambda s: s.replace(",",""))

# MERGE BY TRACT: TractFIPS20 -> LocationID
merged = pd.merge(
    cdc,
    ruca,
    left_on='LocationID',
    right_on='TractFIPS20',
    how='inner'
)

# CAST LOCATIONID TO STR
merged['LocationID'] = merged['LocationID'].astype(str)

# DROP ROWS WITH RUCA==10 (most rural, lowest commute category)
merged = merged[merged['PrimaryRUCA'] != 10]

# DROP NAN ROWS
merged = merged.dropna()

# FILTER TO STATES WE NEED
midwest = {"Missouri","Illinois","Indiana","Iowa","Kansas","Michigan",
	"Minnesota","Nebraska","North Dakota","South Dakota","Ohio","Wisconsin"}
merged = merged[merged['StateDesc'].isin(midwest)]

# RESET INDEX
merged = merged.reset_index(drop=True)

# LOG
print("FINAL DATA SHAPE")
print(merged.shape)

# SAVE TO FILE
merged.to_csv('../shapes/filtered_diabetes.csv', index=False)

