# Medidas 2 MagicScript

## Tests
Test zone for script development, you can see the readme [here](/Tests/README.md). When script has significant changes, they will be added to [MagicScript.m](MagicScript.m)

## Checkpoints

✅ run simulations on OpenEMS with structures made by us. (at [this](https://github.com/monteiroman/medidas2_magicScript/tree/e91d77f7ba519339ee20ab937bb6875e94559fc0) moment)

✅ Obtain similar horn structure as [1].

❌ Obtain similar simulation values as [1].

❌ Simulate a full corrugated horn.

❌ Test other horn shapes (Hyperbolic, Sine-squared,    Exponential, Tangential).

❌ Make functions for different Horn designs in octave.



## OpenEMS installation on Linux
Install openEMS as detailed [here](http://www.openems.de/index.php/Compile_from_Source.html#Linux) for the 
first instalation. Don't forget to install Paraview (again, it is described on the link).

Next instalations/updates can be done as described 
[here](https://github.com/thliebig/openEMS-Project#update-instruction)

## References
[1] **Dual-Polarization and Low-Sidelobe Corrugated Rectangular Horn Antennas for Outdoor RCS Measurement** Changying Wu*, Congxiang Li, Chufeng Hu and Yevhen Yashchyshyn.

[2] **Design of Corrugated Horns: A Primer** Christophe Granet and Graeme L. James.

[3] **Antenna Designer’s Notebook** [(link)](http://antennadesigner.org/)

[4] **Corrugated Horn Antenna Tutorial for openEMS** [(link)](https://openems.de/forum/viewtopic.php?f=3&t=900)

[5] **Corrugated Horn Antenna Design in MATLAB and CST** [(link)](https://www.youtube.com/watch?v=Fh7Ri-CNEjs&ab_channel=SimulationMaster). Or you can see a modified script derived from that example [here](/MatlabToCST_example).
