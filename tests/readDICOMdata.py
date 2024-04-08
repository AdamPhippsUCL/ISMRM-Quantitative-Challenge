# Test script to read in DICOM files

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
import SimpleITK as sitk
from scipy.io import savemat


# Import DICOM script from imgtools
sys.path.insert(0, r"C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python")
import DICOM# type: ignore



# Define list of image names
imagefnames = glob.glob(r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Data\00005-MR-ep2d_basic_clinical_ORIG\*")

# Create DICOM object
testdcm = DICOM.MRIDICOM(DICOM_fnames = imagefnames, multiframe = False)


# Image array
ImageArray = testdcm.constructImageArray()


# Separate by b value
b0slices = np.asarray(range(0,10))
b500slices = np.asarray(range(10,40))
b1000slices = np.asarray(range(40,70))
b2000slices = np.asarray(range(70,100))

b0Image = ImageArray[b0slices]
b0Image = np.vstack((b0Image, b0Image, b0Image))

b500Image = ImageArray[b500slices]
b1000Image = ImageArray[b1000slices]
b2000Image = ImageArray[b2000slices]

# Normalise images

normb500Image = b500Image/b0Image
normb500Image[np.isnan(normb500Image)] = 0
normb500Image[np.isinf(normb500Image)] = 0
print(np.sum(np.isnan(normb500Image)))

normb1000Image = b1000Image/b0Image
normb1000Image[np.isnan(normb1000Image)] = 0
normb1000Image[np.isinf(normb1000Image)] = 0

normb2000Image = b2000Image/b0Image
normb2000Image[np.isnan(normb2000Image)] = 0
normb2000Image[np.isinf(normb2000Image)] = 0

# == Save normalised images as mha

ImageFolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\Images"

sitk.WriteImage(sitk.GetImageFromArray(b0Image), f'{ImageFolder}/b0.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb500Image), f'{ImageFolder}/normb500.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb1000Image), f'{ImageFolder}/normb1000.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb2000Image), f'{ImageFolder}/normb2000.mha')

# plt.figure()
# plt.imshow(normb500Image[5], vmax = 1, cmap = 'gray')
# plt.colorbar()
# plt.show()