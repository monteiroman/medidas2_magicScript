'# MWS Version: Version 2021.1 - Nov 10 2020 - ACIS 30.0.1 -

'# length = mm
'# frequency = GHz
'# time = ns
'# frequency range: fmin = frequency_minimum fmax = frequency_maximum


'@ use template: Antenna (Horn, Waveguide)

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
' Template for Antenna in Free Space
' ==================================
' (CSTxMWSxONLY)
' draw the bounding box
Plot.DrawBox True
' set units to mm, ghz
With Units 
     .Geometry "mm"
     .Frequency "ghz"
     .Time "ns" 
End With 
' set background material to vacuum
With Background 
     .Type "Normal" 
     .Epsilon "1.0" 
     .Mue "1.0" 
     .XminSpace "0.0" 
     .XmaxSpace "0.0" 
     .YminSpace "0.0" 
     .YmaxSpace "0.0" 
     .ZminSpace "0.0" 
     .ZmaxSpace "0.0" 
End With 
' set boundary conditions to open
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "none" 
     .Ysymmetry "none" 
     .Zsymmetry "none" 
End With
' switch on FD-TET setting for accurate farfields
FDSolver.ExtrudeOpenBC "True" 
Mesh.FPBAAvoidNonRegUnite "True" 
Mesh.ConsiderSpaceForLowerMeshLimit "False" 
Mesh.MinimumStepNumber "5"

'@ new component: Antenna

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Component.New "Antenna"

'@ define brick: Antenna:Waveguide_section

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Brick
     .Reset 
     .Name "Waveguide_section" 
     .Component "Antenna" 
     .Material "PEC" 
     .Xrange "-waveguide_height/2", "waveguide_height/2" 
     .Yrange "-waveguide_width/2", "waveguide_width/2" 
     .Zrange "-waveguide_length", "0" 
     .Create
End With

'@ pick face

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.PickFaceFromId "Antenna:Waveguide_section", "2"

'@ define port: 1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Port 
     .Reset 
     .PortNumber "1" 
     .NumberOfModes "1" 
     .AdjustPolarization False 
     .PolarizationAngle "0.0" 
     .ReferencePlaneDistance "0" 
     .TextSize "50" 
     .Coordinates "Picks" 
     .Orientation "positive" 
     .PortOnBound "True" 
     .ClipPickedPortToBound "False" 
     .Xrange "-0.58835", "0.58835" 
     .Yrange "-1.1767", "1.1767" 
     .Zrange "-4.4969", "-4.4969" 
     .XrangeAdd "0.0", "0.0" 
     .YrangeAdd "0.0", "0.0" 
     .ZrangeAdd "0.0", "0.0" 
     .SingleEnded "False" 
     .Create 
End With

'@ pick face

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.PickFaceFromId "Antenna:Waveguide_section", "2"

'@ define material: portbackPEC

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Material 
     .Reset 
     .Name "portbackPEC"
     .FrqType "all" 
     .Type "Pec" 
     .Epsilon "1.0" 
     .Mue "1.0" 
     .Rho "0" 
     .ThermalType "Normal" 
     .ThermalConductivity "0" 
     .HeatCapacity "0" 
     .MetabolicRate "0" 
     .BloodFlow "0" 
     .Colour "0.647059", "0.666667", "0.72549" 
     .Wireframe "True" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .Create
End With

'@ define extrude: Antenna:port-back

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Extrude 
     .Reset 
     .Name "port-back" 
     .Component "Antenna" 
     .Material "portbackPEC" 
     .Mode "Picks" 
     .Height "waveguide_length/20" 
     .Twist "0.0" 
     .Taper "0.0" 
     .UsePicksForHeight "False" 
     .DeleteBaseFaceSolid "False" 
     .ClearPickedFace "True" 
     .Create 
End With

'@ define brick: Antenna:Front_face

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Brick
     .Reset 
     .Name "Front_face" 
     .Component "Antenna" 
     .Material "Vacuum" 
     .Xrange "-aperture_height/2", "aperture_height/2" 
     .Yrange "-aperture_width/2", "aperture_width/2" 
     .Zrange "flare_length", "flare_length" 
     .Create
End With

'@ new curve: Corner_outline

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Curve.NewCurve "Corner_outline"

'@ activate local coordinates

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
WCS.ActivateWCS "local"

'@ pick end point

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.PickEndpointFromId "Antenna:Front_face", "2"

'@ pick end point

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.PickEndpointFromId "Antenna:Waveguide_section", "2"

'@ pick end point

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.PickEndpointFromId "Antenna:Waveguide_section", "5"

'@ align wcs with three points

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
WCS.AlignWCSWithSelected "3Points"

'@ store picked point: 1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.NextPickToDatabase "1" 
Pick.PickEndpointFromId "Antenna:Waveguide_section", "5"

'@ snap point to drawplane

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.SnapLastPointToDrawplane

'@ store picked point: 2

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.NextPickToDatabase "2" 
Pick.PickEndpointFromId "Antenna:Waveguide_section", "2"

'@ snap point to drawplane

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.SnapLastPointToDrawplane

'@ store picked point: 3

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Pick.NextPickToDatabase "3" 
Pick.PickEndpointFromId "Antenna:Front_face", "2"

'@ define curve polygon: Corner_outline:Horn

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Polygon 
     .Reset 
     .Name "Horn" 
     .Curve "Corner_outline" 
     .Point "xp(1)", "yp(1)" 
     .LineTo "xp(2)", "yp(2)" 
     .LineTo "xp(3)", "yp(3)" 
     .Create 
End With

'@ define curve arc: Corner_outline:Rounding

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Arc
     .Reset 
     .Name "Rounding" 
     .Curve "Corner_outline" 
     .Orientation "Clockwise" 
     .XCenter "0" 
     .YCenter "rounding_radius" 
     .X1 "0.0" 
     .Y1 "0.0" 
     .X2 "0.0" 
     .Y2 "0.0" 
     .Angle "arc_sweep_angle" 
     .UseAngle "True" 
     .Segments "0" 
     .Create
End With

'@ activate global coordinates

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
WCS.ActivateWCS "global"

'@ transform curve: mirror Corner_outline

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Horn

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Horn_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "0", "1", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .MirrorCurve 
End With

'@ delete shapes

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.Delete "Antenna:Front_face" 
Solid.Delete "Antenna:Waveguide_section"

'@ define curveloft: Antenna:Rounding_section_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With LoftCurves
     .Reset 
     .Name "Rounding_section_1" 
     .Component "Antenna" 
     .Material "PEC" 
     .Solid "False" 
     .MinimizeTwist "True" 
     .AddCurve "Corner_outline:Horn_1_1" 
     .AddCurve "Corner_outline:Horn_1" 
     .Create 
End With

'@ transform curve: mirror Corner_outline:Horn

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Horn_2

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding_2

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding_2" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .MirrorCurve 
End With

'@ define curveloft: Antenna:Rounding_section_2

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With LoftCurves
     .Reset 
     .Name "Rounding_section_2" 
     .Component "Antenna" 
     .Material "PEC" 
     .Solid "False" 
     .MinimizeTwist "True" 
     .AddCurve "Corner_outline:Rounding" 
     .AddCurve "Corner_outline:Rounding_2" 
     .Create 
End With

'@ transform curve: mirror Corner_outline:Horn_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Horn_2_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Horn_2_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "True" 
     .Destination "" 
     .MirrorCurve 
End With

'@ transform curve: mirror Corner_outline:Rounding_2_1

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Transform 
     .Reset 
     .Name "Corner_outline:Rounding_2_1" 
     .Origin "Free" 
     .Center "0", "0", "0" 
     .PlaneNormal "1", "0", "0" 
     .MultipleObjects "True" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Destination "" 
     .MirrorCurve 
End With

'@ define curveloft: Antenna:Rounding_section_3

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With LoftCurves
     .Reset 
     .Name "Rounding_section_3" 
     .Component "Antenna" 
     .Material "PEC" 
     .Solid "False" 
     .MinimizeTwist "True" 
     .AddCurve "Corner_outline:Horn_1" 
     .AddCurve "Corner_outline:Horn_1_1" 
     .Create 
End With

'@ define curveloft: Antenna:Rounding_section_4

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With LoftCurves
     .Reset 
     .Name "Rounding_section_4" 
     .Component "Antenna" 
     .Material "PEC" 
     .Solid "False" 
     .MinimizeTwist "True" 
     .AddCurve "Corner_outline:Horn_2_1_1" 
     .AddCurve "Corner_outline:Horn_2_1" 
     .Create 
End With

'@ boolean add shapes: Antenna:Rounding_section_1, Antenna:Rounding_section_2

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.Add "Antenna:Rounding_section_1", "Antenna:Rounding_section_2"

'@ boolean add shapes: Antenna:Rounding_section_1, Antenna:Rounding_section_3

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.Add "Antenna:Rounding_section_1", "Antenna:Rounding_section_3"

'@ boolean add shapes: Antenna:Rounding_section_1, Antenna:Rounding_section_4

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.Add "Antenna:Rounding_section_1", "Antenna:Rounding_section_4"

'@ rename block: Antenna:Rounding_section_1 to: Antenna:Complete_antenna

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.Rename "Antenna:Rounding_section_1", "Complete_antenna"

'@ thicken sheet: Antenna:Complete_antenna

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
If (metal_thickness > 0) Then Solid.ThickenSheetAdvanced "Antenna:Complete_antenna", "Outside", "metal_thickness", "True"

'@ define frequency range

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solver.FrequencyRange "frequency_minimum", "frequency_maximum"

'@ define boundaries

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "electric" 
     .Ysymmetry "magnetic" 
     .Zsymmetry "none" 
     .XminThermal "isothermal" 
     .XmaxThermal "isothermal" 
     .YminThermal "isothermal" 
     .YmaxThermal "isothermal" 
     .ZminThermal "isothermal" 
     .ZmaxThermal "isothermal" 
     .XsymmetryThermal "isothermal" 
     .YsymmetryThermal "isothermal" 
     .ZsymmetryThermal "none" 
     .ApplyInAllDirections "False" 
     .XminTemperature "" 
     .XminTemperatureType "None" 
     .XmaxTemperature "" 
     .XmaxTemperatureType "None" 
     .YminTemperature "" 
     .YminTemperatureType "None" 
     .YmaxTemperature "" 
     .YmaxTemperatureType "None" 
     .ZminTemperature "" 
     .ZminTemperatureType "None" 
     .ZmaxTemperature "" 
     .ZmaxTemperatureType "None" 
End With

'@ define solver parameters

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .CalculationType "TD-S" 
     .StimulationPort "All" 
     .StimulationMode "All" 
     .SteadyStateLimit "-40" 
     .MeshAdaption "False" 
     .AutoNormImpedance "False" 
     .NormingImpedance "50" 
     .CalculateModesOnly "False" 
     .SParaSymmetry "False" 
     .StoreTDResultsInCache "False" 
     .FullDeembedding "False" 
     .UseDistributedComputing "False" 
     .DistributeMatrixCalculation "False" 
     .MPIParallelization "False" 
     .SuperimposePLWExcitation "False" 
End With

'@ execute macro: Multiple_farfield_requests

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Dim f_min_norm As Double, f_max_norm As Double
Dim cnt As Integer
Dim f_monitor_norm As Double
Dim monitor_name As String
MakeSureParameterExists("num_ff_monitors", "0")
'Normalised frequency range
f_min_norm = frequency_minimum/frequency_centre
f_max_norm = frequency_maximum/frequency_centre
For cnt=1 To num_ff_monitors
	f_monitor_norm = f_min_norm + (cnt-1)*(f_max_norm-f_min_norm)/(num_ff_monitors-1)
                monitor_name ="ff_"+Replace(Format(f_monitor_norm,"0.000"),",",".")+"xf0"
	With Monitor
	     .Reset
	     .Name monitor_name
	     .FieldType "Farfield"
	     .Frequency Replace(CStr(f_monitor_norm), ",", ".") +"*frequency_centre"
	     .Create
	End With
Next cnt

'@ define automesh parameters

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
With Mesh 
     .AutomeshStraightLines "True" 
     .AutomeshEllipticalLines "True" 
     .AutomeshRefineAtPecLines "True", "2" 
     .AutomeshRefinePecAlongAxesOnly "False" 
     .AutomeshAtEllipseBounds "True", "10" 
     .AutomeshAtWireEndPoints "True" 
     .AutomeshAtProbePoints "True" 
     .SetAutomeshRefineDielectricsType "Generalized" 
     .MergeThinPECLayerFixpoints "True" 
     .EquilibrateMesh "False" 
     .EquilibrateMeshRatio "1.19" 
     .UseCellAspectRatio "False" 
     .CellAspectRatio "50.0" 
     .UsePecEdgeModel "True" 
     .MeshType "PBA" 
     .AutoMeshLimitShapeFaces "True" 
     .AutoMeshNumberOfShapeFaces "1000" 
     .PointAccEnhancement "0" 
     .SurfaceOptimization "True" 
     .SurfaceSmoothing "3" 
     .MinimumCurvatureRefinement "100" 
     .CurvatureRefinementFactor "0.05" 
     .AnisotropicCurvatureRefinement "False" 
     .SmallFeatureSize "0.0" 
     .SurfaceTolerance "0.0" 
     .SurfaceToleranceType "Relative" 
     .NormalTolerance "22.5" 
     .AnisotropicCurvatureRefinementFSM "False" 
     .SurfaceMeshEnrichment "0" 
     .DensityTransitionsFSM "0.5" 
     .VolumeOptimization "True" 
     .VolumeSmoothing "True" 
     .VolumeMeshMethod "Delaunay" 
     .SurfaceMeshMethod "General" 
     .SurfaceMeshGeometryAccuracy "1.0e-6" 
     .DelaunayOptimizationLevel "2" 
     .DelaunayPropagationFactor "1.050000" 
     .DensityTransitions "0.5" 
     .MeshAllRegions "False" 
     .ConvertGeometryDataAfterMeshing "True" 
     .AutomeshFixpointsForBackground "True" 
     .PBAType "Fast PBA" 
     .AutomaticPBAType "True" 
     .FPBAAvoidNonRegUnite "True" 
     .DetectSmallSolidPEC "False" 
     .ConsiderSpaceForLowerMeshLimit "True" 
     .RatioLimitGovernsLocalRefinement "False" 
     .GapDetection "False" 
     .FPBAGapTolerance "1e-3" 
     .MaxParallelThreads "8"
     .SetParallelMode "Maximum"
End With 
With Solver 
     .UseSplitComponents "True" 
     .PBAFillLimit "75" 
     .EnableSubgridding "False" 
     .AlwaysExcludePec "False" 
End With

'@ deactivate fast model update

'[VERSION]2009.8|18.0.3|20090230[/VERSION]
Solid.FastModelUpdate "False"

'@ set units in materials

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
Material.SetUnitInMaterial "$CoilMaterial$", "GHz", "mm" 
Material.SetUnitInMaterial "portbackPEC", "GHz", "mm"

'@ set mesh properties (for backward compatibility)

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
With MeshSettings
     .SetMeshType "Hex"
     .Set "Version", 0%
     .SetMeshType "Tet"
     .Set "Version", 0%
     .SetMeshType "Srf"
     .Set "Version", 0%
End With
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "CellsPerWavelengthPolicy", "cellsperwavelength" 
     .Set "CurvatureOrderPolicy", "off" 
     .SetMeshType "Plane" 
     .Set "CurvatureOrderPolicy", "off" 
End With

'@ change solver type

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
ChangeSolverType "HF Time Domain"

'@ set mesh properties (Hexahedral)

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
With Mesh 
     .MeshType "PBA" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "10" 
     .Set "StepsPerWaveFar", "10" 
     .Set "WavelengthRefinementSameAsNear", "1" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "10" 
     .Set "StepsPerBoxFar", "1" 
     .Set "MaxStepNear", "0" 
     .Set "MaxStepFar", "0" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     .Set "UseMaxStepAbsolute", "0" 
     .Set "GeometryRefinementSameAsNear", "0" 
     'MIN CELL 
     .Set "UseRatioLimitGeometry", "1" 
     .Set "RatioLimitGeometry", "10" 
     .Set "MinStepGeometryX", "0" 
     .Set "MinStepGeometryY", "0" 
     .Set "MinStepGeometryZ", "0" 
     .Set "UseSameMinStepGeometryXYZ", "1" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "FaceRefinementOn", "0" 
     .Set "FaceRefinementPolicy", "2" 
     .Set "FaceRefinementRatio", "2" 
     .Set "FaceRefinementStep", "0" 
     .Set "FaceRefinementNSteps", "2" 
     .Set "EllipseRefinementOn", "0" 
     .Set "EllipseRefinementPolicy", "2" 
     .Set "EllipseRefinementRatio", "2" 
     .Set "EllipseRefinementStep", "0" 
     .Set "EllipseRefinementNSteps", "2" 
     .Set "FaceRefinementBufferLines", "3" 
     .Set "EdgeRefinementOn", "1" 
     .Set "EdgeRefinementPolicy", "1" 
     .Set "EdgeRefinementRatio", "2" 
     .Set "EdgeRefinementStep", "0" 
     .Set "EdgeRefinementBufferLines", "3" 
     .Set "RefineEdgeMaterialGlobal", "0" 
     .Set "RefineAxialEdgeGlobal", "0" 
     .Set "BufferLinesNear", "3" 
     .Set "UseDielectrics", "1" 
     .Set "EquilibrateOn", "0" 
     .Set "Equilibrate", "1.5" 
     .Set "IgnoreThinPanelMaterial", "0" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "SnapToAxialEdges", "1"
     .Set "SnapToPlanes", "1"
     .Set "SnapToSpheres", "1"
     .Set "SnapToEllipses", "1"
     .Set "SnapToCylinders", "1"
     .Set "SnapToCylinderCenters", "1"
     .Set "SnapToEllipseCenters", "1"
End With 
With Discretizer 
     .MeshType "PBA" 
     .PBAType "Fast PBA" 
     .AutomaticPBAType "True" 
     .FPBAAccuracyEnhancement "enable"
     .ConnectivityCheck "False"
     .ConvertGeometryDataAfterMeshing "True" 
     .UsePecEdgeModel "True" 
     .GapDetection "False" 
     .FPBAGapTolerance "1e-3" 
     .SetMaxParallelMesherThreads "Hex", "8"
     .SetParallelMesherMode "Hex", "Maximum"
     .PointAccEnhancement "0" 
     .UseSplitComponents "True" 
     .EnableSubgridding "False" 
     .PBAFillLimit "75" 
     .AlwaysExcludePec "False" 
End With

'@ define time domain solver parameters

'[VERSION]2014.6|23.0.0|20141128[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .Method "Hexahedral"
     .CalculationType "TD-S"
     .StimulationPort "All"
     .StimulationMode "All"
     .SteadyStateLimit "-40"
     .MeshAdaption "False"
     .AutoNormImpedance "False"
     .NormingImpedance "50"
     .CalculateModesOnly "False"
     .SParaSymmetry "False"
     .StoreTDResultsInCache  "False"
     .FullDeembedding "False"
     .SuperimposePLWExcitation "False"
     .UseSensitivityAnalysis "False"
End With

'@ ___________Magus Anchor Point Creation___________

'[VERSION]2014.6|23.0.0|20090230[/VERSION]

            
        '[VERSION]2017.5|26.0.1|20170804[/VERSION]


        '@ execute macro: CreateMagusBoundaryAnchorPoints

        '[VERSION]2017.5|26.0.1|20170804[/VERSION]
        ' Extract and store the original boundary conditions
        Dim XminBoundaryType As String
        Dim XmaxBoundaryType As String
        Dim YminBoundaryType As String
        Dim YmaxBoundaryType As String
        Dim ZminBoundaryType As String
        Dim ZmaxBoundaryType As String
        XminBoundaryType = Boundary.GetXmin
        XmaxBoundaryType = Boundary.GetXmax
        YminBoundaryType = Boundary.GetYmin
        YmaxBoundaryType = Boundary.GetYmax
        ZminBoundaryType = Boundary.GetZmin
        ZmaxBoundaryType = Boundary.GetZmax
        ' Set the boundary to open in all directions
        Dim TightBoundary As String
        TightBoundary = "open"
        Boundary.Xmin(TightBoundary)
        Boundary.Xmax(TightBoundary)
        Boundary.Ymin(TightBoundary)
        Boundary.Ymax(TightBoundary)
        Boundary.Zmin(TightBoundary)
        Boundary.Zmax(TightBoundary)
        ' Setup WCS
        WCS.AlignWCSWithGlobalCoordinates
        ' Variables used to extract the "tight" bounding box
        Dim xmin As Double
        Dim xmax As Double
        Dim ymin As Double
        Dim ymax As Double
        Dim zmin As Double
        Dim zmax As Double
        Boundary.GetCalculationBox(xmin, xmax, ymin, ymax, zmin, zmax)
        ' Set the the original boundary conditions back to what they were
        Boundary.Xmin(XminBoundaryType)
        Boundary.Xmax(XmaxBoundaryType)
        Boundary.Ymin(YminBoundaryType)
        Boundary.Ymax(YmaxBoundaryType)
        Boundary.Zmin(ZminBoundaryType)
        Boundary.Zmax(ZmaxBoundaryType)
        ' Anchor Points
        Dim AnchorPointFolderName As String
        AnchorPointFolderName = "Boundary"
        ' Move the local WCS to the boundaries, orientate so that the normal points outward
        Dim du As Double
        Dim dv As Double
        Dim dw As Double
        ' Zmax
        du = 0.5*( xmin + xmax)
        dv = 0.5*( ymin + ymax)
        dw = zmax
        WCS.MoveWCS "local", du, dv, dw
        ' Set Zmax Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Zmax"
        ' Store WCS of current Anchor Point
        WCS.Store "Zmax"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        ' Zmin
        du = 0.5*( xmin + xmax)
        dv = 0.5*( ymin + ymax)
        dw = zmin
        WCS.MoveWCS "local", du, dv, dw
        WCS.SetNormal "0", "0", "-1"
        ' Set Zmin Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Zmin"
        ' Store WCS of current Anchor Point
        WCS.Store "Zmin"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        ' Xmax
        du = xmax
        dv = 0.5*( ymin + ymax)
        dw = 0.5*( zmin + zmax)
        WCS.MoveWCS "local", du, dv, dw
        WCS.SetNormal "1", "0", "0"
        ' Set Zmin Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Xmax"
        ' Store WCS of current Anchor Point
        WCS.Store "Xmax"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        ' Xmin
        du = xmin
        dv = 0.5*( ymin + ymax)
        dw = 0.5*( zmin + zmax)
        WCS.MoveWCS "local", du, dv, dw
        WCS.SetNormal "-1", "0", "0"
        ' Set Zmin Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Xmin"
        ' Store WCS of current Anchor Point
        WCS.Store "Xmin"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        ' Ymax
        du = 0.5*( xmin + xmax)
        dv = ymax
        dw = 0.5*( zmin + zmax)
        WCS.MoveWCS "local", du, dv, dw
        WCS.SetNormal "0", "1", "0"
        ' Set Zmin Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Ymax"
        ' Store WCS of current Anchor Point
        WCS.Store "Ymax"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        ' Ymin
        du = 0.5*( xmin + xmax)
        dv = ymin
        dw = 0.5*( zmin + zmax)
        WCS.MoveWCS "local", du, dv, dw
        WCS.SetNormal "0", "-1", "0"
        ' Set Zmin Anchor Point
        AnchorPoint.Store AnchorPointFolderName + ":" + "Ymin"
        ' Store WCS of current Anchor Point
        WCS.Store "Ymin"
        'Realign WCS with global
        WCS.AlignWCSWithGlobalCoordinates
        
      '@ ___________End: Magus Anchor Point Creation___________

'@ Exported from Antenna Magus: Aperture-matched waveguide-fed pyramidal horn - Friday, July 9, 2021

'[VERSION]2014.6|23.0.0|20090230[/VERSION]


'@ set mesh properties (for backward compatibility)

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With MeshSettings 
    .SetMeshType "Hex"
    .Set "PlaneMergeVersion", 0
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "8" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "8" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "8" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ set mesh properties (Hexahedral)

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With Mesh 
     .MeshType "PBA" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "5" 
     .Set "StepsPerWaveFar", "5" 
     .Set "WavelengthRefinementSameAsNear", "1" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "5" 
     .Set "StepsPerBoxFar", "1" 
     .Set "MaxStepNear", "0" 
     .Set "MaxStepFar", "0" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     .Set "UseMaxStepAbsolute", "0" 
     .Set "GeometryRefinementSameAsNear", "0" 
     'MIN CELL 
     .Set "UseRatioLimitGeometry", "1" 
     .Set "RatioLimitGeometry", "10" 
     .Set "MinStepGeometryX", "0" 
     .Set "MinStepGeometryY", "0" 
     .Set "MinStepGeometryZ", "0" 
     .Set "UseSameMinStepGeometryXYZ", "1" 
End With 
With MeshSettings 
     .Set "PlaneMergeVersion", "2" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "FaceRefinementOn", "0" 
     .Set "FaceRefinementPolicy", "2" 
     .Set "FaceRefinementRatio", "2" 
     .Set "FaceRefinementStep", "0" 
     .Set "FaceRefinementNSteps", "2" 
     .Set "EllipseRefinementOn", "0" 
     .Set "EllipseRefinementPolicy", "2" 
     .Set "EllipseRefinementRatio", "2" 
     .Set "EllipseRefinementStep", "0" 
     .Set "EllipseRefinementNSteps", "2" 
     .Set "FaceRefinementBufferLines", "3" 
     .Set "EdgeRefinementOn", "1" 
     .Set "EdgeRefinementPolicy", "1" 
     .Set "EdgeRefinementRatio", "2" 
     .Set "EdgeRefinementStep", "0" 
     .Set "EdgeRefinementBufferLines", "3" 
     .Set "RefineEdgeMaterialGlobal", "0" 
     .Set "RefineAxialEdgeGlobal", "0" 
     .Set "BufferLinesNear", "3" 
     .Set "UseDielectrics", "1" 
     .Set "EquilibrateOn", "0" 
     .Set "Equilibrate", "1.5" 
     .Set "IgnoreThinPanelMaterial", "0" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "SnapToAxialEdges", "1"
     .Set "SnapToPlanes", "1"
     .Set "SnapToSpheres", "1"
     .Set "SnapToEllipses", "1"
     .Set "SnapToCylinders", "1"
     .Set "SnapToCylinderCenters", "1"
     .Set "SnapToEllipseCenters", "1"
End With 
With Discretizer 
     .ConnectivityCheck "False"
     .UsePecEdgeModel "True" 
     .PointAccEnhancement "0" 
     .TSTVersion "0"
	  .PBAVersion "2009023009" 
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "8" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "False" 
     .ShowStructureProfile "False" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ set mesh properties (Hexahedral)

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With Mesh 
     .MeshType "PBA" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "15" 
     .Set "StepsPerWaveFar", "15" 
     .Set "WavelengthRefinementSameAsNear", "1" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "15" 
     .Set "StepsPerBoxFar", "1" 
     .Set "MaxStepNear", "0" 
     .Set "MaxStepFar", "0" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     .Set "UseMaxStepAbsolute", "0" 
     .Set "GeometryRefinementSameAsNear", "0" 
     'MIN CELL 
     .Set "UseRatioLimitGeometry", "1" 
     .Set "RatioLimitGeometry", "10" 
     .Set "MinStepGeometryX", "0" 
     .Set "MinStepGeometryY", "0" 
     .Set "MinStepGeometryZ", "0" 
     .Set "UseSameMinStepGeometryXYZ", "1" 
End With 
With MeshSettings 
     .Set "PlaneMergeVersion", "2" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "FaceRefinementOn", "0" 
     .Set "FaceRefinementPolicy", "2" 
     .Set "FaceRefinementRatio", "2" 
     .Set "FaceRefinementStep", "0" 
     .Set "FaceRefinementNSteps", "2" 
     .Set "EllipseRefinementOn", "0" 
     .Set "EllipseRefinementPolicy", "2" 
     .Set "EllipseRefinementRatio", "2" 
     .Set "EllipseRefinementStep", "0" 
     .Set "EllipseRefinementNSteps", "2" 
     .Set "FaceRefinementBufferLines", "3" 
     .Set "EdgeRefinementOn", "1" 
     .Set "EdgeRefinementPolicy", "1" 
     .Set "EdgeRefinementRatio", "2" 
     .Set "EdgeRefinementStep", "0" 
     .Set "EdgeRefinementBufferLines", "3" 
     .Set "RefineEdgeMaterialGlobal", "0" 
     .Set "RefineAxialEdgeGlobal", "0" 
     .Set "BufferLinesNear", "3" 
     .Set "UseDielectrics", "1" 
     .Set "EquilibrateOn", "0" 
     .Set "Equilibrate", "1.5" 
     .Set "IgnoreThinPanelMaterial", "0" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "SnapToAxialEdges", "1"
     .Set "SnapToPlanes", "1"
     .Set "SnapToSpheres", "1"
     .Set "SnapToEllipses", "1"
     .Set "SnapToCylinders", "1"
     .Set "SnapToCylinderCenters", "1"
     .Set "SnapToEllipseCenters", "1"
End With 
With Discretizer 
     .ConnectivityCheck "False"
     .UsePecEdgeModel "True" 
     .PointAccEnhancement "0" 
     .TSTVersion "0"
	  .PBAVersion "2009023009" 
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "8" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "True" 
     .ShowStructureProfile "True" 
     .SetStructureTransparent "True" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ set mesh properties (Hexahedral)

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With Mesh 
     .MeshType "PBA" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "5" 
     .Set "StepsPerWaveFar", "5" 
     .Set "WavelengthRefinementSameAsNear", "1" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "5" 
     .Set "StepsPerBoxFar", "1" 
     .Set "MaxStepNear", "0" 
     .Set "MaxStepFar", "0" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     .Set "UseMaxStepAbsolute", "0" 
     .Set "GeometryRefinementSameAsNear", "0" 
     'MIN CELL 
     .Set "UseRatioLimitGeometry", "1" 
     .Set "RatioLimitGeometry", "10" 
     .Set "MinStepGeometryX", "0" 
     .Set "MinStepGeometryY", "0" 
     .Set "MinStepGeometryZ", "0" 
     .Set "UseSameMinStepGeometryXYZ", "1" 
End With 
With MeshSettings 
     .Set "PlaneMergeVersion", "2" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "FaceRefinementOn", "0" 
     .Set "FaceRefinementPolicy", "2" 
     .Set "FaceRefinementRatio", "2" 
     .Set "FaceRefinementStep", "0" 
     .Set "FaceRefinementNSteps", "2" 
     .Set "EllipseRefinementOn", "0" 
     .Set "EllipseRefinementPolicy", "2" 
     .Set "EllipseRefinementRatio", "2" 
     .Set "EllipseRefinementStep", "0" 
     .Set "EllipseRefinementNSteps", "2" 
     .Set "FaceRefinementBufferLines", "3" 
     .Set "EdgeRefinementOn", "1" 
     .Set "EdgeRefinementPolicy", "1" 
     .Set "EdgeRefinementRatio", "2" 
     .Set "EdgeRefinementStep", "0" 
     .Set "EdgeRefinementBufferLines", "3" 
     .Set "RefineEdgeMaterialGlobal", "0" 
     .Set "RefineAxialEdgeGlobal", "0" 
     .Set "BufferLinesNear", "3" 
     .Set "UseDielectrics", "1" 
     .Set "EquilibrateOn", "0" 
     .Set "Equilibrate", "1.5" 
     .Set "IgnoreThinPanelMaterial", "0" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "SnapToAxialEdges", "1"
     .Set "SnapToPlanes", "1"
     .Set "SnapToSpheres", "1"
     .Set "SnapToEllipses", "1"
     .Set "SnapToCylinders", "1"
     .Set "SnapToCylinderCenters", "1"
     .Set "SnapToEllipseCenters", "1"
End With 
With Discretizer 
     .ConnectivityCheck "False"
     .UsePecEdgeModel "True" 
     .PointAccEnhancement "0" 
     .TSTVersion "0"
	  .PBAVersion "2009023009" 
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "1" 
     .Step2 "1" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "True" 
     .ShowStructureProfile "True" 
     .SetStructureTransparent "True" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

'@ set mesh properties (Hexahedral)

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With Mesh 
     .MeshType "PBA" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "10" 
     .Set "StepsPerWaveFar", "10" 
     .Set "WavelengthRefinementSameAsNear", "1" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "10" 
     .Set "StepsPerBoxFar", "1" 
     .Set "MaxStepNear", "0" 
     .Set "MaxStepFar", "0" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     .Set "UseMaxStepAbsolute", "0" 
     .Set "GeometryRefinementSameAsNear", "0" 
     'MIN CELL 
     .Set "UseRatioLimitGeometry", "1" 
     .Set "RatioLimitGeometry", "10" 
     .Set "MinStepGeometryX", "0" 
     .Set "MinStepGeometryY", "0" 
     .Set "MinStepGeometryZ", "0" 
     .Set "UseSameMinStepGeometryXYZ", "1" 
End With 
With MeshSettings 
     .Set "PlaneMergeVersion", "2" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "FaceRefinementOn", "0" 
     .Set "FaceRefinementPolicy", "2" 
     .Set "FaceRefinementRatio", "2" 
     .Set "FaceRefinementStep", "0" 
     .Set "FaceRefinementNSteps", "2" 
     .Set "EllipseRefinementOn", "0" 
     .Set "EllipseRefinementPolicy", "2" 
     .Set "EllipseRefinementRatio", "2" 
     .Set "EllipseRefinementStep", "0" 
     .Set "EllipseRefinementNSteps", "2" 
     .Set "FaceRefinementBufferLines", "3" 
     .Set "EdgeRefinementOn", "1" 
     .Set "EdgeRefinementPolicy", "1" 
     .Set "EdgeRefinementRatio", "2" 
     .Set "EdgeRefinementStep", "0" 
     .Set "EdgeRefinementBufferLines", "3" 
     .Set "RefineEdgeMaterialGlobal", "0" 
     .Set "RefineAxialEdgeGlobal", "0" 
     .Set "BufferLinesNear", "3" 
     .Set "UseDielectrics", "1" 
     .Set "EquilibrateOn", "0" 
     .Set "Equilibrate", "1.5" 
     .Set "IgnoreThinPanelMaterial", "0" 
End With 
With MeshSettings 
     .SetMeshType "Hex" 
     .Set "SnapToAxialEdges", "1"
     .Set "SnapToPlanes", "1"
     .Set "SnapToSpheres", "1"
     .Set "SnapToEllipses", "1"
     .Set "SnapToCylinders", "1"
     .Set "SnapToCylinderCenters", "1"
     .Set "SnapToEllipseCenters", "1"
End With 
With Discretizer 
     .ConnectivityCheck "False"
     .UsePecEdgeModel "True" 
     .PointAccEnhancement "0" 
     .TSTVersion "0"
	  .PBAVersion "2009023009" 
End With

'@ farfield plot options

'[VERSION]2021.1|30.0.1|20201110[/VERSION]
With FarfieldPlot 
     .Plottype "Cartesian" 
     .Vary "angle1" 
     .Theta "0" 
     .Phi "0" 
     .Step "5" 
     .Step2 "5" 
     .SetLockSteps "True" 
     .SetPlotRangeOnly "False" 
     .SetThetaStart "0" 
     .SetThetaEnd "180" 
     .SetPhiStart "0" 
     .SetPhiEnd "360" 
     .SetTheta360 "True" 
     .SymmetricRange "True" 
     .SetTimeDomainFF "False" 
     .SetFrequency "5" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "True" 
     .ShowStructureProfile "True" 
     .SetStructureTransparent "True" 
     .SetFarfieldTransparent "False" 
     .AspectRatio "Free" 
     .ShowGridlines "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "Gain" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .IncludeUnitCellSidewalls "True" 
     .SetScaleLinear "False" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .SetMaxReferenceMode "abs" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .SetAntennaType "unknown" 
     .Phistart "1.000000e+00", "0.000000e+00", "0.000000e+00" 
     .Thetastart "0.000000e+00", "0.000000e+00", "1.000000e+00" 
     .PolarizationVector "0.000000e+00", "1.000000e+00", "0.000000e+00" 
     .SetCoordinateSystemType "spherical" 
     .SetAutomaticCoordinateSystem "True" 
     .SetPolarizationType "Linear" 
     .SlantAngle 0.000000e+00 
     .Origin "bbox" 
     .Userorigin "0.000000e+00", "0.000000e+00", "0.000000e+00" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+00" 
     .LossyGround "False" 
     .GroundEpsilon "1" 
     .GroundKappa "0" 
     .EnablePhaseCenterCalculation "False" 
     .SetPhaseCenterAngularLimit "3.000000e+01" 
     .SetPhaseCenterComponent "boresight" 
     .SetPhaseCenterPlane "both" 
     .ShowPhaseCenter "True" 
     .ClearCuts 

     .StoreSettings
End With

