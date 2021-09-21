# Testing structures
***
## [Structure_test](https://github.com/monteiroman/medidas2_magicScript/tree/main/Tests/Structure_test)

### 1) test.m
Simple extruded horn face. It was a test in how to make a straight horn wall from the calculated corrugated profile.

### 2) test_walls.m
Semi-corrugated horn inside a mesh.

### 3) test_horn_with_chokes.m
Simulation of horn designed on [1] (see references in this [readme](/README.md)).

<img src="ReadmeData/test_horn_with_chokes/openEMS_Structure_chokes.png" alt="openEMS_Structure_chokes.png" width="800"/>

***
## [Simulation_test](https://github.com/monteiroman/medidas2_magicScript/tree/main/Tests/Simulation_test)

### 1) test_first.m
First test of a piramidal corrugated horn.

### 2) test_second.m
Test with piramidal semi-corrugated horn. At the moment we are facing some memory problems with this script but some results can be seen below.

|2D Structure parts obtained from Octave|3D Structure in paraView|
|:---:|:---:|
|<img src="ReadmeData/test_second/openEMS_Structure.png" alt="openEMS_Structure.png" width="800"/>|<img src="ReadmeData/test_second/paraView_Structure.png" alt="paraView_Structure.png" width="800"/>|

|Reflection coeficient|Farfield Directivity|
|:---:|:---:|
|<img src="ReadmeData/test_second/openEMS_Reflection_coeficient.png" alt="openEMS_Reflection_coeficient.png" width="800"/>|<img src="ReadmeData/test_second/openEMS_Farfield_Directivity.png" alt="openEMS_Farfield_Directivity.png" width="800"/>|

|Farfield Directivity polar coordinates|Farfield Directivity Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/test_second/openEMS_Farfield_Directivity_Polar.png" alt="openEMS_Farfield_Directivity_Polar.png" width="800"/>|<img src="ReadmeData/test_second/openEMS_Farfield_Directivity_Ludwig3.png" alt="openEMS_Farfield_Directivity_Ludwig3.png" width="800"/>|

|Radiation Patern OpenEMS|Radiation Patern and 3D structure in Paraview|
|:---:|:---:|
|<img src="ReadmeData/test_second/openEMS_Radiation_patern.png" alt="openEMS_Radiation_patern.png" width="800"/>|<img src="ReadmeData/test_second/paraView_Radiation_patern.png" alt="paraView_Radiation_patern.png" width="800"/>|


### 3) test_horn_with_chokes_simulation.m

Test with corrugated horn designed at [1](see references [here](/README.md)). At the moment we are facing some memory problems with this script but some results can be seen below.

|2D Structure parts obtained from Octave|3D Structure in paraView|
|:---:|:---:|
|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Structure.png" alt="openEMS_Structure.png" width="800"/>|<img src="ReadmeData/test_horn_with_chokes_simulation/paraView_Structure.png" alt="paraView_Structure.png" width="800"/>|

|Reflection coeficient|Farfield Directivity|
|:---:|:---:|
|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Reflection_coeficient.png" alt="openEMS_Reflection_coeficient.png" width="800"/>|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Farfield_Directivity.png" alt="openEMS_Farfield_Directivity.png" width="800"/>|

|Farfield Directivity polar coordinates|Farfield Directivity Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Farfield_Directivity_Polar.png" alt="openEMS_Farfield_Directivity_Polar.png" width="800"/>|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Farfield_Directivity_Ludwig3.png" alt="openEMS_Farfield_Directivity_Ludwig3.png" width="800"/>|

|Radiation Patern OpenEMS|Radiation Patern and 3D structure in Paraview|
|:---:|:---:|
|<img src="ReadmeData/test_horn_with_chokes_simulation/openEMS_Radiation_patern.png" alt="openEMS_Radiation_patern.png" width="800"/>|<img src="ReadmeData/test_horn_with_chokes_simulation/paraView_Radiation_patern.png" alt="paraView_Radiation_patern.png" width="800"/>|

### 4) test_corrugated_horn.m

Tangential profile Corrugated Horn

|2D Structure parts obtained from Octave|3D Structure in paraView|
|:---:|:---:|
|<img src="ReadmeData/test_corrugated_horn/openEMS_Structure_tan.png" alt="openEMS_Structure_tan.png" width="800"/>|<img src="ReadmeData/test_corrugated_horn/paraView_Structure_tan.png" alt="paraView_Structure_tan.png" width="800"/>|

|Reflection coeficient|Farfield Directivity|
|:---:|:---:|
|<img src="ReadmeData/test_corrugated_horn/openEMS_Reflection_coeficient_tan.png" alt="openEMS_Reflection_coeficient_tan.png" width="800"/>|<img src="ReadmeData/test_corrugated_horn/openEMS_Farfield_Directivity_tan.png" alt="openEMS_Farfield_Directivity_tan.png" width="800"/>|

|Farfield Directivity polar coordinates|Farfield Directivity Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/test_corrugated_horn/openEMS_Farfield_Directivity_Polar_tan.png" alt="openEMS_Farfield_Directivity_Polar_tan.png" width="800"/>|<img src="ReadmeData/test_corrugated_horn/openEMS_Farfield_Directivity_Ludwig3_tan.png" alt="openEMS_Farfield_Directivity_Ludwig3_tan.png" width="800"/>|

|Radiation Patern OpenEMS|Radiation Patern and 3D structure in Paraview|
|:---:|:---:|
|<img src="ReadmeData/test_corrugated_horn/openEMS_Radiation_patern_tan.png" alt="openEMS_Radiation_patern_tan.png" width="800"/>|<img src="ReadmeData/test_corrugated_horn/paraView_Radiation_patern_tan.png" alt="paraView_Radiation_patern_tan.png" width="800"/>|