import sys
import config as cfg
import data

## Inputs:
# PathToData: Training Data
# Num: Number of training samples
# StartNum: Start processing from this sample number
## Outputs:
# voxels_preprocessed.vtu

dataPath = cfg.Data_path
num = cfg.num
startNum = cfg.startNum

valid = 0
for i in range( startNum, num+startNum ):
    print(str(i) + "/" + str(num))
    if data.preprocess( dataPath, i ):
        valid += 1

print( "Converted {:d} samples.".format( valid ) )
