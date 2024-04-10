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

b500xslices = np.asarray(range(10,20))
b500yslices = np.asarray(range(20,30))
b500zslices = np.asarray(range(30,40))

b1000xslices = np.asarray(range(40,50))
b1000yslices = np.asarray(range(50,60))
b1000zslices = np.asarray(range(60,70))

b2000xslices = np.asarray(range(70,80))
b2000yslices = np.asarray(range(80,90))
b2000zslices = np.asarray(range(90,100))

b0Image = ImageArray[b0slices]

b500xImage = ImageArray[b500xslices]
b500yImage = ImageArray[b500yslices]
b500zImage = ImageArray[b500zslices]


b1000xImage = ImageArray[b1000xslices]
b1000yImage = ImageArray[b1000yslices]
b1000zImage = ImageArray[b1000zslices]

b2000xImage = ImageArray[b2000xslices]
b2000yImage = ImageArray[b2000yslices]
b2000zImage = ImageArray[b2000zslices]


# Normalise images

normb500xImage = b500xImage/b0Image
normb500xImage[np.isnan(normb500xImage)] = 0
normb500xImage[np.isinf(normb500xImage)] = 0
normb500yImage = b500yImage/b0Image
normb500yImage[np.isnan(normb500yImage)] = 0
normb500yImage[np.isinf(normb500yImage)] = 0
normb500zImage = b500zImage/b0Image
normb500zImage[np.isnan(normb500zImage)] = 0
normb500zImage[np.isinf(normb500zImage)] = 0

normb1000xImage = b1000xImage/b0Image
normb1000xImage[np.isnan(normb1000xImage)] = 0
normb1000xImage[np.isinf(normb1000xImage)] = 0
normb1000yImage = b1000yImage/b0Image
normb1000yImage[np.isnan(normb1000yImage)] = 0
normb1000yImage[np.isinf(normb1000yImage)] = 0
normb1000zImage = b1000zImage/b0Image
normb1000zImage[np.isnan(normb1000zImage)] = 0
normb1000zImage[np.isinf(normb1000zImage)] = 0

normb2000xImage = b2000xImage/b0Image
normb2000xImage[np.isnan(normb2000xImage)] = 0
normb2000xImage[np.isinf(normb2000xImage)] = 0
normb2000yImage = b2000yImage/b0Image
normb2000yImage[np.isnan(normb2000yImage)] = 0
normb2000yImage[np.isinf(normb2000yImage)] = 0
normb2000zImage = b2000zImage/b0Image
normb2000zImage[np.isnan(normb2000zImage)] = 0
normb2000zImage[np.isinf(normb2000zImage)] = 0

# == Save normalised images as mha

ImageFolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\Images"

sitk.WriteImage(sitk.GetImageFromArray(b0Image), f'{ImageFolder}/b0.mha')

sitk.WriteImage(sitk.GetImageFromArray(normb500xImage), f'{ImageFolder}/normb500x.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb500yImage), f'{ImageFolder}/normb500y.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb500zImage), f'{ImageFolder}/normb500z.mha')

sitk.WriteImage(sitk.GetImageFromArray(normb1000xImage), f'{ImageFolder}/normb1000x.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb1000yImage), f'{ImageFolder}/normb1000y.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb1000zImage), f'{ImageFolder}/normb1000z.mha')

sitk.WriteImage(sitk.GetImageFromArray(normb2000xImage), f'{ImageFolder}/normb2000x.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb2000yImage), f'{ImageFolder}/normb2000y.mha')
sitk.WriteImage(sitk.GetImageFromArray(normb2000zImage), f'{ImageFolder}/normb2000z.mha')

# plt.figure()
# plt.imshow(normb500Image[5], vmax = 1, cmap = 'gray')
# plt.colorbar()
# plt.show()