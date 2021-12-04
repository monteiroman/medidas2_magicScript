# Medidas 2 MagicScript
This project was developed for the Department of Medidas Electrónicas II (Electronic Measurements II) of the UTN FRBA (National Technological University of Buenos Aires).

## MagicScript
The script can design and simulate rectangular horn antennas defined by the user. It also can sweep horn parameters automatically and stores results on "outputs/".

Parameters that can be swept are:

||Parameter|Description|
|:---:|:---:|:---:|
|1|ai|Input width|
|2|bi|Input height|
|3|ao|Output width|
|4|bo|Output height|
|5|corr_step|Corrugation step in mm|
|6|delta|Ratio between corrugation slot and teeth|
|7|depth_a|A wall corrugation depth|
|8|depth_b|B wall corrugation depth|
|9|wg_length|Input waveguide length|
|10|num_of_corrugations|Number of horn flare corrugations|
|11|a_jump|First jump in A wall corrugations|
|12|b_jump|First jump in B wall corrugations|

⚠️**NOTE**: This script only runs on Linux. OpenEMS must be installed on "~/opt/openEMS" (as it is by default) in order to use CSXCad.

***
## Results

**Obtained results with OpenEMS (no "b jump")**

Results shown below are for [these](/ReadmeData/no_b_jump/0_horn_variables.txt) horn dimensions.

|3D Structure in Paraview|3D Radiation patern in Paraview|
|:---:|:---:|
|<img src="ReadmeData/no_b_jump/paraview_Structure.png" alt="paraview_Structure.png" width="800"/>|<img src="ReadmeData/no_b_jump/paraview_Farfield.png" alt="paraview_Farfield.png" width="800"/>|

|Reflection coeficient|Farfield polar coordinates|
|:---:|:---:|
|<img src="ReadmeData/no_b_jump/0_S11_output.png" alt="0_S11_output.png" width="800"/>|<img src="ReadmeData/no_b_jump/0_Farfield_directivity_polar_output.png" alt="0_Farfield_directivity_polar_output.png" width="800"/>|

|Farfield|Farfield Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/no_b_jump/0_Farfield_directivity_output.png" alt="0_Farfield_directivity_output.png" width="800"/>|<img src="ReadmeData/no_b_jump/0_Farfield_directivity_ludwig_output.png" alt="0_Farfield_directivity_ludwig_output.png" width="800"/>|

|Horn Structure|
|:---:|
|<img src="ReadmeData/no_b_jump/0_structure_output.png" alt="0_structure_output.png" width="600"/>|

**Obtained results with OpenEMS (1mm "b jump")**

Results shown below are for [these](/ReadmeData/1mm_b_jump/1_horn_variables.txt) horn dimensions.

|3D Structure in Paraview|3D Radiation patern in Paraview|
|:---:|:---:|
|<img src="ReadmeData/1mm_b_jump/paraview_Structure.png" alt="paraview_Structure.png" width="800"/>|<img src="ReadmeData/1mm_b_jump/paraview_Farfield.png" alt="paraview_Farfield.png" width="800"/>|

|Reflection coeficient|Farfield polar coordinates|
|:---:|:---:|
|<img src="ReadmeData/1mm_b_jump/1_S11_output.png" alt="1_S11_output.png" width="800"/>|<img src="ReadmeData/1mm_b_jump/1_Farfield_directivity_polar_output.png" alt="1_Farfield_directivity_polar_output.png" width="800"/>|

|Farfield|Farfield Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/1mm_b_jump/1_Farfield_directivity_output.png" alt="1_Farfield_directivity_output.png" width="800"/>|<img src="ReadmeData/1mm_b_jump/1_Farfield_directivity_ludwig_output.png" alt="1_Farfield_directivity_ludwig_output.png" width="800"/>|

|Horn Structure|
|:---:|
|<img src="ReadmeData/1mm_b_jump/1_structure_output.png" alt="1_structure_output.png" width="600"/>|

**Obtained results with OpenEMS (5mm "b jump")**

Results shown below are for [these](/ReadmeData/5mm_b_jump/9_horn_variables.txt) horn dimensions.

|3D Structure in Paraview|3D Radiation patern in Paraview|
|:---:|:---:|
|<img src="ReadmeData/5mm_b_jump/paraview_Structure.png" alt="paraview_Structure.png" width="800"/>|<img src="ReadmeData/5mm_b_jump/paraview_Farfield.png" alt="paraview_Farfield.png" width="800"/>|

|Reflection coeficient|Farfield polar coordinates|
|:---:|:---:|
|<img src="ReadmeData/5mm_b_jump/9_S11_output.png" alt="9_S11_output.png" width="800"/>|<img src="ReadmeData/5mm_b_jump/9_Farfield_directivity_polar_output.png" alt="9_Farfield_directivity_polar_output.png" width="800"/>|

|Farfield|Farfield Ludwig3 coordinates|
|:---:|:---:|
|<img src="ReadmeData/5mm_b_jump/9_Farfield_directivity_output.png" alt="9_Farfield_directivity_output.png" width="800"/>|<img src="ReadmeData/5mm_b_jump/9_Farfield_directivity_ludwig_output.png" alt="9_Farfield_directivity_ludwig_output.png" width="800"/>|

|Horn Structure|
|:---:|
|<img src="ReadmeData/5mm_b_jump/9_structure_output.png" alt="9_structure_output.png" width="600"/>|

***
## Project State
|Step|Description|Status|
|:---:|:---:|:---:|
|1|Run simulations on OpenEMS with structures made by us.|✅|
|2|Obtain similar horn structure as [1].|✅|
|3|Obtain similar simulation values as [1].|✅|
|4|Simulate a full corrugated horn.|✅|
|5|Test other horn shapes (Hyperbolic, Exponential, Tangential).|✅|
|6|Make function for Horn design in octave.|✅|
|7|Make a parameter sweep script for horn antennas in octave.|✅|
|8|Waveguide to N connector adapter|✅|

***
## Tests
Project tests can be seen [here](Tests/), they were the first steps developing the root project. You can see the readme [here](/Tests/README.md).

***
## OpenEMS installation on Linux
Install openEMS as detailed [here](http://www.openems.de/index.php/Compile_from_Source.html#Linux) for the 
first instalation. Don't forget to install Paraview (again, it is described on the link).

You may face some problems with 'libvtk5-dev' and 'libvtk5-qt4-dev'. Please refer to [this](https://ubuntuforums.org/showthread.php?t=2431395) post to solve them.

Next instalations/updates can be done as described 
[here](https://github.com/thliebig/openEMS-Project#update-instruction).

***
## References
[1] **Dual-Polarization and Low-Sidelobe Corrugated Rectangular Horn Antennas for Outdoor RCS Measurement** Changying Wu, Congxiang Li, Chufeng Hu and Yevhen Yashchyshyn.

[2] **Design of Corrugated Horns: A Primer** Christophe Granet and Graeme L. James.

[3] **Antenna Designer’s Notebook** [(link)](http://antennadesigner.org/).

[4] **Corrugated Horn Antenna Tutorial for openEMS** [(link)](https://openems.de/forum/viewtopic.php?f=3&t=900).

[5] **Corrugated Horn Antenna Design in MATLAB and CST** [(link)](https://www.youtube.com/watch?v=Fh7Ri-CNEjs&ab_channel=SimulationMaster). Or you can see a modified script derived from that example [here](/MatlabToCST_example).

[6] **Antenna Theory: Analysis and Design** Constantine A. Balanis
