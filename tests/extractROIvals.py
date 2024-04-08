# Import relevant libraries

import numpy as np
import matplotlib.pyplot as plt
import os
import SimpleITK as sitk
from scipy.io import savemat

# Define image folder
ImageFolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\Images"

# Define ROI folder
ROIFolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROIs"


# Define b value
b=2000

# Define ROI number
ROINum = 5


# Read image
bImage = sitk.GetArrayFromImage(sitk.ReadImage(f'{ImageFolder}/normb{b}.mha'))
# bImage = sitk.GetArrayFromImage(sitk.ReadImage(f'{ImageFolder}/b{b}.mha'))

# Read ROI
ROI = sitk.GetArrayFromImage(sitk.ReadImage(f'{ROIFolder}/ROI{ROINum}.img'))
ROI = (ROI==1)

# Extract ROI values
ROIvals = bImage[ROI]


# == Save ROI vals as mat

ROIValsFolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROI Values"


try:
    os.makedirs(f'{ROIValsFolder}/b{b}')
except:
    None

savemat(
    f'{ROIValsFolder}/b{b}/ROI{ROINum}.mat',
    dict(zip(['ROIvals'], [ROIvals]))
)

plt.figure()
plt.hist(ROIvals)
plt.show()