The bash script nondicom_writer reads a DICOM file of a 4d MRI scan, extracts velocity and mag information, and then writes that information in a concise and convenient file format (together with some relevant metadata) to the directory nondicom_data. This script calls Matlab functions in the subdirectory MatlabScripts

The folder nondicom_reader contains c++ code to reads the info written by nondicom_writer and represents it in a C++ "grid" object defined in the core tools directory. The grid object has (will soon have) functionality to export to .vtk files; a format that can be displayed with VisIt and other software. 


