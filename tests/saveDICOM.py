import pydicom as pyd


DICOMfname = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Data\00005-MR-ep2d_basic_clinical_ORIG\00042.dcm"


dcm = pyd.dcmread(DICOMfname)


with open( f'dcm.txt', 'w') as f:
    f.write(str(dcm))