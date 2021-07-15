'# MWS Version: Version 2015.0 - Jan 16 2015 - ACIS 24.0.2 -

'# length = in
'# frequency = GHz
'# time = ns
'# frequency range: fmin = 8 fmax = 12


'@ use template: Antenna (in Free Space, waveguide)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
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

'@ define units

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Units 
     .Geometry "in" 
     .Frequency "GHz" 
     .Time "ns" 
     .Temperature "Celsius" 
     .Voltage "V" 
     .Current "A" 
     .Resistance "Ohm" 
     .Conductance "S" 
     .Capacitance "pF" 
     .Inductance "nH" 
End With

'@ new component: component1

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Component.New "component1"

'@ define brick: component1:solid1

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Brick
     .Reset 
     .Name "solid1" 
     .Component "component1" 
     .Material "PEC" 
     .Xrange "-0", "1" 
     .Yrange "-0", "0.5" 
     .Zrange "0", "0.5" 
     .Create
End With

'@ pick face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickFaceFromId "component1:solid1", "1"

'@ align wcs with face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
WCS.AlignWCSWithSelectedFace 
Pick.PickCenterpointFromId "component1:solid1", "1" 
WCS.AlignWCSWithSelectedPoint

'@ move wcs

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
WCS.MoveWCS "local", "0.0", "0.0", "2"

'@ define cylinder: component1:solid2

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Cylinder 
     .Reset 
     .Name "solid2" 
     .Component "component1" 
     .Material "PEC" 
     .OuterRadius "r1" 
     .InnerRadius "0" 
     .Axis "z" 
     .Zrange "0", "0.25" 
     .Xcenter "0" 
     .Ycenter "0" 
     .Segments "0" 
     .Create 
End With

'@ switch working plane

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Plot.DrawWorkplane "false"

'@ pick face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickFaceFromId "component1:solid2", "1"

'@ pick face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickFaceFromId "component1:solid1", "1"

'@ define loft: component1:solid3

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Loft 
     .Reset 
     .Name "solid3" 
     .Component "component1" 
     .Material "PEC" 
     .Tangency "0.150000" 
     .CreateNew 
End With

'@ clear picks

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.ClearAllPicks

'@ boolean add shapes: component1:solid1, component1:solid2

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Solid.Add "component1:solid1", "component1:solid2"

'@ boolean add shapes: component1:solid1, component1:solid3

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Solid.Add "component1:solid1", "component1:solid3"

'@ pick face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickFaceFromId "component1:solid1", "9"

'@ pick face

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickFaceFromId "component1:solid1", "11"

'@ shell object: component1:solid1

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Solid.AdvancedShell "component1:solid1", "Outside", "0.02"

'@ pick end point

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickEndpointFromId "component1:solid1", "15"

'@ pick edge

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Pick.PickEdgeFromId "component1:solid1", "21", "17"

'@ define port: 1

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
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
     .Xrange "0", "1" 
     .Yrange "0", "0.5" 
     .Zrange "0", "0" 
     .XrangeAdd "0.0", "0.0" 
     .YrangeAdd "0.0", "0.0" 
     .ZrangeAdd "0.0", "0.0" 
     .SingleEnded "False" 
     .Create 
End With

'@ activate global coordinates

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
WCS.ActivateWCS "global"

'@ define frequency range

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Solver.FrequencyRange "8", "12"

'@ define boundaries

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Boundary
     .Xmin "expanded open" 
     .Xmax "expanded open" 
     .Ymin "expanded open" 
     .Ymax "expanded open" 
     .Zmin "expanded open" 
     .Zmax "expanded open" 
     .Xsymmetry "magnetic" 
     .Ysymmetry "electric" 
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

'@ define monitor: e-field (f=10)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Monitor 
     .Reset 
     .Name "e-field (f=10)" 
     .Dimension "Volume" 
     .Domain "Frequency" 
     .FieldType "Efield" 
     .Frequency "10" 
     .Create 
End With

'@ define farfield monitor: farfield (f=10)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Monitor 
     .Reset 
     .Name "farfield (f=10)" 
     .Domain "Frequency" 
     .FieldType "Farfield" 
     .Frequency "10" 
     .Create 
End With

'@ set mesh properties

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Mesh 
     .UseRatioLimit "True" 
     .RatioLimit "5" 
     .LinesPerWavelength "5" 
     .MinimumStepNumber "5" 
     .Automesh "True" 
     .MeshType "PBA" 
End With

'@ set shape accuracy

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Solid.ShapeVisualizationAccuracy "96" 
Solid.ShapeVisualizationOffset "25" 
Pick.ClearAllPicks

'@ define probe: E-field (Farfield) (Phi; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Probe
     .Reset 
     .Name "E-field (Farfield) (Phi; 0 0 200)" 
     .Field "EFarfield" 
     .Orientation "Phi" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define probe: H-field (Farfield) (Theta; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Probe
     .Reset 
     .Name "H-field (Farfield) (Theta; 0 0 200)" 
     .Field "HFarfield" 
     .Orientation "Theta" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define probe: E-field (Farfield) (Theta; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Probe.Delete "E-field (Farfield) (Phi; 0 0 200)" 
With Probe
     .Reset 
     .Name "E-field (Farfield) (Theta; 0 0 200)" 
     .Field "EFarfield" 
     .Orientation "Theta" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define probe: H-field (Farfield) (Phi; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Probe.Delete "H-field (Farfield) (Theta; 0 0 200)" 
With Probe
     .Reset 
     .Name "H-field (Farfield) (Phi; 0 0 200)" 
     .Field "HFarfield" 
     .Orientation "Phi" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define probe: E-field (Farfield) (Phi; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Probe.Delete "E-field (Farfield) (Theta; 0 0 200)" 
With Probe
     .Reset 
     .Name "E-field (Farfield) (Phi; 0 0 200)" 
     .Field "EFarfield" 
     .Orientation "Phi" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define probe: H-field (Farfield) (Theta; 0 0 200)

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
Probe.Delete "H-field (Farfield) (Phi; 0 0 200)" 
With Probe
     .Reset 
     .Name "H-field (Farfield) (Theta; 0 0 200)" 
     .Field "HFarfield" 
     .Orientation "Theta" 
     .SetPosition1 "0" 
     .SetPosition2 "0" 
     .SetPosition3 "200" 
     .SetSampling "16.0" 
     .SetCoordinateSystemType "Spherical" 
     .Create
End With

'@ define automesh parameters

'[VERSION]2010.0|20.0.0|20090230[/VERSION]
With Mesh 
     .AutomeshStraightLines "True" 
     .AutomeshEllipticalLines "True" 
     .AutomeshRefineAtPecLines "True", "2" 
     .AutomeshRefinePecAlongAxesOnly "False" 
     .AutomeshAtEllipseBounds "True", "10" 
     .AutomeshAtWireEndPoints "True" 
     .AutomeshAtProbePoints "True" 
     .SetAutomeshRefineDielectricsType "Generalized" 
     .MergeThinPECLayerFixpoints "False" 
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
     .AutomaticPBAType "False" 
     .FPBAAvoidNonRegUnite "True" 
     .DetectSmallSolidPEC "False" 
     .ConsiderSpaceForLowerMeshLimit "True" 
     .RatioLimitGovernsLocalRefinement "False" 
     .GapDetection "False" 
     .FPBAGapTolerance "1e-3" 
     .SetMaxParallelMesherThreads "Hex", "8"
     .SetParallelMesherMode "Hex", "Maximum"
     .AutomeshRefineThermalMaterials "False" 
     .SetThermalRefinementConductivityReference "1e-3" 
     .SetThermalRefinementHeatCapacityReference "1e-3" 
     .SetParallelMesherMode "Tet", "maximum" 
     .SetMaxParallelMesherThreads "Tet", "1" 
     .ConnectivityCheck "False"
End With 
With Solver 
     .UseSplitComponents "True" 
     .PBAFillLimit "99" 
     .EnableSubgridding "False" 
     .AlwaysExcludePec "False" 
End With

'@ set mesh properties (Hexahedral)

'[VERSION]2014.0|23.0.0|20130901[/VERSION]
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
     .Set "RefineEdgesAtBoundary", "0" 
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
     .AutomaticPBAType "False" 
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
     .PBAFillLimit "99" 
     .AlwaysExcludePec "False" 
End With

'@ define solver parameters

'[VERSION]2011.0|21.0.0|20100920[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .CalculationType "TD-S"
     .StimulationPort "All"
     .StimulationMode "All"
     .SteadyStateLimit "-30"
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

'@ define time domain solver parameters

'[VERSION]2012.0|22.0.0|20111003[/VERSION]
Mesh.SetCreator "High Frequency" 
With Solver 
     .Method "Hexahedral"
     .CalculationType "TD-S"
     .StimulationPort "All"
     .StimulationMode "All"
     .SteadyStateLimit "-30"
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

'@ farfield plot options

'[VERSION]2012.0|22.0.0|20111003[/VERSION]
With FarfieldPlot 
     .Plottype "3D" 
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
     .SetFrequency "10" 
     .SetTime "0" 
     .SetColorByValue "True" 
     .DrawStepLines "False" 
     .DrawIsoLongitudeLatitudeLines "False" 
     .ShowStructure "True" 
     .SetStructureTransparent "False" 
     .SetFarfieldTransparent "True" 
     .SetSpecials "enablepolarextralines" 
     .SetPlotMode "EField" 
     .Distance "1" 
     .UseFarfieldApproximation "True" 
     .SetScaleLinear "True" 
     .SetLogRange "40" 
     .SetLogNorm "0" 
     .DBUnit "0" 
     .EnableFixPlotMaximum "False" 
     .SetFixPlotMaximumValue "1.0" 
     .SetInverseAxialRatio "False" 
     .SetAxesType "user" 
     .Phistart "1.000000e+000", "0.000000e+000", "0.000000e+000" 
     .Thetastart "0.000000e+000", "0.000000e+000", "1.000000e+000" 
     .PolarizationVector "0.000000e+000", "1.000000e+000", "0.000000e+000" 
     .SetCoordinateSystemType "ludwig3" 
     .SetPolarizationType "Linear" 
     .SlantAngle 4.500000e+001 
     .Origin "free" 
     .Userorigin "5.000000e-001", "2.500000e-001", "2.602000e+000" 
     .SetUserDecouplingPlane "False" 
     .UseDecouplingPlane "False" 
     .DecouplingPlaneAxis "X" 
     .DecouplingPlanePosition "0.000000e+000" 
     .EnablePhaseCenterCalculation "True" 
     .SetPhaseCenterAngularLimit "1.500000e+001" 
     .SetPhaseCenterComponent "phi" 
     .SetPhaseCenterPlane "H-plane" 
     .ShowPhaseCenter "True" 
     .StoreSettings
End With

'@ set mesh properties

'[VERSION]2013.0|23.0.0|20130115[/VERSION]
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "CellsPerWavelengthPolicy", "cellsperwavelength" 
     .Set "CurvatureOrderPolicy", "off" 
     .SetMeshType "Plane" 
     .Set "CurvatureOrderPolicy", "off" 
End With

'@ change solver type

'[VERSION]2013.0|23.0.0|20130115[/VERSION]
ChangeSolverType "HF Time Domain"

'@ set pba mesh type

'[VERSION]2014.0|23.0.0|20130901[/VERSION]
Mesh.MeshType "PBA"

'@ define pml specials

'[VERSION]2014.0|23.0.0|20130901[/VERSION]
With Boundary
     .ReflectionLevel "0.0001" 
     .MinimumDistanceType "Fraction" 
     .MinimumDistancePerWavelengthNewMeshEngine "4" 
     .MinimumDistanceReferenceFrequencyType "Center" 
     .FrequencyForMinimumDistance "10" 
     .SetAbsoluteDistance "0.0" 
End With

'@ set pba mesh type

'[VERSION]2014.0|23.0.0|20130901[/VERSION]
Mesh.MeshType "PBA"

'@ set mesh properties (for backward compatibility)

'[VERSION]2015.0|24.0.2|20140727[/VERSION]
With MeshSettings 
    .SetMeshType "Hex"
    .Set "PlaneMergeVersion", 0
End With

'@ define background

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Background 
     .ResetBackground 
     .XminSpace "10" 
     .XmaxSpace "10" 
     .YminSpace "10" 
     .YmaxSpace "10" 
     .ZminSpace "10" 
     .ZmaxSpace "10" 
     .ApplyInAllDirections "True" 
End With 
With Material 
     .Reset 
     .FrqType "all"
     .Type "Normal"
     .MaterialUnit "Frequency", "Hz"
     .MaterialUnit "Geometry", "m"
     .MaterialUnit "Time", "s"
     .MaterialUnit "Temperature", "Kelvin"
     .Epsilon "1.0"
     .Mue "1.0"
     .Sigma "0.0"
     .TanD "0.0"
     .TanDFreq "0.0"
     .TanDGiven "False"
     .TanDModel "ConstSigma"
     .EnableUserConstTanDModelOrderEps "False"
     .ConstTanDModelOrderEps "1"
     .SetElParametricConductivity "False"
     .ReferenceCoordSystem "Global"
     .CoordSystemType "Cartesian"
     .SigmaM "0"
     .TanDM "0.0"
     .TanDMFreq "0.0"
     .TanDMGiven "False"
     .TanDMModel "ConstSigma"
     .EnableUserConstTanDModelOrderMue "False"
     .ConstTanDModelOrderMue "1"
     .SetMagParametricConductivity "False"
     .DispModelEps  "None"
     .DispModelMue "None"
     .DispersiveFittingSchemeEps "Nth Order"
     .MaximalOrderNthModelFitEps "10"
     .ErrorLimitNthModelFitEps "0.1"
     .UseOnlyDataInSimFreqRangeNthModelEps "False"
     .DispersiveFittingSchemeMue "Nth Order"
     .MaximalOrderNthModelFitMue "10"
     .ErrorLimitNthModelFitMue "0.1"
     .UseOnlyDataInSimFreqRangeNthModelMue "False"
     .UseGeneralDispersionEps "False"
     .UseGeneralDispersionMue "False"
     .NLAnisotropy "False"
     .NLAStackingFactor "1"
     .NLADirectionX "1"
     .NLADirectionY "0"
     .NLADirectionZ "0"
     .Rho "0.0"
     .ThermalType "Normal"
     .ThermalConductivity "0.0"
     .HeatCapacity "0.0"
     .MetabolicRate "0"
     .BloodFlow "0"
     .VoxelConvection "0"
     .MechanicsType "Unused"
     .Colour "0.6", "0.6", "0.6" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .ChangeBackgroundMaterial
End With

'@ define material: Iron

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Material
     .Reset
     .Name "Iron"
     .Folder ""
.FrqType "hf" 
.Type "Lossy metal" 
.SetMaterialUnit "Hz", "mm" 
.Mue "1" 
.Kappa "1.04e+007" 
.FrqType "all" 
.Type "Nonlinear" 
.SetMaterialUnit "Hz", "mm" 
.Mue "1000" 
.NonlinearMeasurementError "1e-3" 
.Kappa "1.04e+007"
.ResetHBList
.AddHBValue " 0.000000e+000", " 0.000000e+000" 
.AddHBValue " 1.592000e+001", " 4.500000e-002" 
.AddHBValue " 3.183000e+001", " 4.900000e-001" 
.AddHBValue " 4.775000e+001", " 7.800000e-001" 
.AddHBValue " 6.366000e+001", " 9.900000e-001" 
.AddHBValue " 7.958000e+001", " 1.130000e+000" 
.AddHBValue " 1.592000e+002", " 1.420000e+000" 
.AddHBValue " 3.183000e+002", " 1.570000e+000" 
.AddHBValue " 4.775000e+002", " 1.621000e+000" 
.AddHBValue " 6.366000e+002", " 1.641000e+000" 
.AddHBValue " 7.958000e+002", " 1.656000e+000" 
.AddHBValue " 1.592000e+003", " 1.697000e+000" 
.AddHBValue " 3.183000e+003", " 1.764000e+000" 
.AddHBValue " 4.775000e+003", " 1.801000e+000" 
.AddHBValue " 6.366000e+003", " 1.838000e+000" 
.AddHBValue " 7.958000e+003", " 1.870000e+000" 
.AddHBValue " 1.592000e+004", " 2.000000e+000" 
.AddHBValue " 3.183000e+004", " 2.136000e+000" 
.AddHBValue " 4.775000e+004", " 2.185000e+000" 
.AddHBValue " 6.366000e+004", " 2.220000e+000" 
.AddHBValue " 7.958000e+004", " 2.250000e+000" 
.AddHBValue " 1.592000e+005", " 2.355000e+000" 
.AddHBValue " 3.183000e+005", " 2.556000e+000" 
.Rho "7870.0" 
.ThermalType "Normal" 
.ThermalConductivity "79.5" 
.HeatCapacity "0.45" 
.MetabolicRate "0" 
.BloodFlow "0" 
.VoxelConvection "0" 
.MechanicsType "Isotropic" 
.YoungsModulus "200" 
.PoissonsRatio "0.291" 
.ThermalExpansionRate "12" 
.FrqType "static" 
.Type "Nonlinear" 
.SetMaterialUnit "Hz", "mm" 
.Mue "1000" 
.NonlinearMeasurementError "1e-3" 
.Kappa "1.04e+007"
.ResetHBList
.AddHBValue " 0.000000e+000", " 0.000000e+000" 
.AddHBValue " 1.592000e+001", " 4.500000e-002" 
.AddHBValue " 3.183000e+001", " 4.900000e-001" 
.AddHBValue " 4.775000e+001", " 7.800000e-001" 
.AddHBValue " 6.366000e+001", " 9.900000e-001" 
.AddHBValue " 7.958000e+001", " 1.130000e+000" 
.AddHBValue " 1.592000e+002", " 1.420000e+000" 
.AddHBValue " 3.183000e+002", " 1.570000e+000" 
.AddHBValue " 4.775000e+002", " 1.621000e+000" 
.AddHBValue " 6.366000e+002", " 1.641000e+000" 
.AddHBValue " 7.958000e+002", " 1.656000e+000" 
.AddHBValue " 1.592000e+003", " 1.697000e+000" 
.AddHBValue " 3.183000e+003", " 1.764000e+000" 
.AddHBValue " 4.775000e+003", " 1.801000e+000" 
.AddHBValue " 6.366000e+003", " 1.838000e+000" 
.AddHBValue " 7.958000e+003", " 1.870000e+000" 
.AddHBValue " 1.592000e+004", " 2.000000e+000" 
.AddHBValue " 3.183000e+004", " 2.136000e+000" 
.AddHBValue " 4.775000e+004", " 2.185000e+000" 
.AddHBValue " 6.366000e+004", " 2.220000e+000" 
.AddHBValue " 7.958000e+004", " 2.250000e+000" 
.AddHBValue " 1.592000e+005", " 2.355000e+000" 
.AddHBValue " 3.183000e+005", " 2.556000e+000" 
.Colour "1", "0", "0" 
.Wireframe "False" 
.Reflection "False" 
.Allowoutline "True" 
.Transparentoutline "False" 
.Transparency "0" 
.Create
End With

'@ define brick: component1:solid2

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Brick
     .Reset 
     .Name "solid2" 
     .Component "component1" 
     .Material "Iron" 
     .Xrange "10", "-10" 
     .Yrange "10", "-10" 
     .Zrange "10", "-10" 
     .Create
End With

'@ define brick: component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Brick
     .Reset 
     .Name "solid3" 
     .Component "component1" 
     .Material "Iron" 
     .Xrange "9.99", "-9.99" 
     .Yrange "9.99", "-9.99" 
     .Zrange "9.99", "-9.99" 
     .Create
End With

'@ boolean subtract shapes: component1:solid2, component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Solid.Subtract "component1:solid2", "component1:solid3"

'@ transform: translate component1:solid2

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "component1:solid2" 
     .Vector "0", "0", "20" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Shape", "Translate" 
End With

'@ pick edge

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Pick.PickEdgeFromId "component1:solid1", "20", "16"

'@ define distance dimension by picks

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Dimension
    .Reset
    .UsePicks True
    .SetType "Distance"
    .SetID "0"
    .SetOrientation "Smart Mode"
    .SetDistance "0.156148"
    .SetViewVector "-0.046720", "-0.060547", "0.997071"
    .Create
End With
Pick.ClearAllPicks

'@ transform: translate component1:solid1

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "component1:solid1" 
     .Vector "-0.5", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Shape", "Translate" 
End With

'@ transform port: translate port1

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "port1" 
     .Vector "-0.5", "0", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Port", "Translate" 
End With

'@ transform: translate component1:solid1

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "component1:solid1" 
     .Vector "0", "-0.25", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Shape", "Translate" 
End With

'@ transform port: translate port1

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "port1" 
     .Vector "0", "-0.25", "0" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Port", "Translate" 
End With

'@ define brick: component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Brick
     .Reset 
     .Name "solid3" 
     .Component "component1" 
     .Material "Iron" 
     .Xrange "-5", "5" 
     .Yrange "-5", "5" 
     .Zrange "10.01", "9.99" 
     .Create
End With

'@ boolean subtract shapes: component1:solid2, component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Solid.Subtract "component1:solid2", "component1:solid3"

'@ define brick: component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Brick
     .Reset 
     .Name "solid3" 
     .Component "component1" 
     .Material "Iron" 
     .Xrange "-5.05", "5.05" 
     .Yrange "-5.05", "5.05" 
     .Zrange "9.99", "10" 
     .Create
End With

'@ transform: translate component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Transform 
     .Reset 
     .Name "component1:solid3" 
     .Vector "0", "0", "-0.010" 
     .UsePickedPoints "False" 
     .InvertPickedPoints "False" 
     .MultipleObjects "False" 
     .GroupObjects "False" 
     .Repetitions "1" 
     .MultipleSelection "False" 
     .Transform "Shape", "Translate" 
End With

'@ delete shape: component1:solid2

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Solid.Delete "component1:solid2"

'@ delete shape: component1:solid3

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Solid.Delete "component1:solid3"

'@ delete dimension 0

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Dimension
    .RemoveDimension "0"
End With

'@ change solver type

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
ChangeSolverType "HF Frequency Domain"

'@ define frequency domain solver parameters

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
Mesh.SetCreator "High Frequency" 
With FDSolver
     .Reset 
     .SetMethod "Tetrahedral", "General purpose" 
     .OrderTet "Second" 
     .OrderSrf "First" 
     .Stimulation "All", "All" 
     .ResetExcitationList 
     .AutoNormImpedance "False" 
     .NormingImpedance "50" 
     .ModesOnly "False" 
     .ConsiderPortLossesTet "True" 
     .SetShieldAllPorts "False" 
     .AccuracyHex "1e-6" 
     .AccuracyTet "1e-4" 
     .AccuracySrf "1e-3" 
     .LimitIterations "False" 
     .MaxIterations "0" 
     .SetCalculateExcitationsInParallel "True", "False", "" 
     .StoreAllResults "False" 
     .StoreResultsInCache "False" 
     .UseHelmholtzEquation "True" 
     .LowFrequencyStabilization "True" 
     .Type "Auto" 
     .MeshAdaptionHex "False" 
     .MeshAdaptionTet "True" 
     .AcceleratedRestart "True" 
     .FreqDistAdaptMode "Distributed" 
     .NewIterativeSolver "True" 
     .TDCompatibleMaterials "False" 
     .ExtrudeOpenBC "True" 
     .SetOpenBCTypeHex "Default" 
     .SetOpenBCTypeTet "Default" 
     .AddMonitorSamples "True" 
     .CalcStatBField "False" 
     .CalcPowerLoss "True" 
     .CalcPowerLossPerComponent "False" 
     .StoreSolutionCoefficients "True" 
     .UseDoublePrecision "False" 
     .UseDoublePrecision_ML "True" 
     .MixedOrderSrf "False" 
     .MixedOrderTet "False" 
     .PreconditionerAccuracyIntEq "0.15" 
     .MLFMMAccuracy "Default" 
     .MinMLFMMBoxSize "0.20" 
     .UseCFIEForCPECIntEq "true" 
     .UseFastRCSSweepIntEq "false" 
     .UseSensitivityAnalysis "False" 
     .SweepErrorThreshold "True", "0.01" 
     .SweepErrorChecks "2" 
     .SweepMinimumSamples "3" 
     .SweepConsiderAll "True" 
     .SweepConsiderReset 
     .SetNumberOfResultDataSamples "1001" 
     .SweepWeightEvanescent "1.0" 
     .AccuracyROM "1e-4" 
     .AddSampleInterval "", "", "1", "Automatic", "True" 
     .AddSampleInterval "", "", "", "Automatic", "False" 
     .MPIParallelization "False"
     .UseDistributedComputing "False"
     .NetworkComputingStrategy "RunRemote"
     .NetworkComputingJobCount "3"
     .LimitCPUs "True"
     .MaxCPUs "48"
End With
With IESolver
     .Reset 
     .UseFastFrequencySweep "True" 
     .UseIEGroundPlane "False" 
     .SetRealGroundMaterialName "" 
     .PreconditionerType "Auto" 
End With
With IESolver
     .SetFMMFFCalcStopLevel "0" 
     .SetFMMFFCalcNumInterpPoints "6" 
     .UseFMMFarfieldCalc "True" 
     .SetCFIEAlpha "0.500000" 
     .LowFrequencyStabilization "False" 
     .LowFrequencyStabilizationML "True" 
     .Multilayer "False" 
     .SetiMoMACC_I "0.0001" 
     .SetiMoMACC_M "0.0001" 
     .DeembedExternalPorts "True" 
     .SetOpenBC_XY "True" 
     .OldRCSSweepDefintion "False" 
     .SetAccuracySetting "Custom" 
     .CalculateSParaforFieldsources "True" 
End With

'@ set mesh properties (Tetrahedral)

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Mesh 
     .MeshType "Tetrahedral" 
     .SetCreator "High Frequency"
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "Version", 1%
     'MAX CELL - WAVELENGTH REFINEMENT 
     .Set "StepsPerWaveNear", "3" 
     .Set "StepsPerWaveFar", "3" 
     .Set "PhaseErrorNear", "0.02" 
     .Set "PhaseErrorFar", "0.02" 
     .Set "CellsPerWavelengthPolicy", "cellsperwavelength" 
     'MAX CELL - GEOMETRY REFINEMENT 
     .Set "StepsPerBoxNear", "5" 
     .Set "StepsPerBoxFar", "1" 
     .Set "ModelBoxDescrNear", "maxedge" 
     .Set "ModelBoxDescrFar", "maxedge" 
     'MIN CELL 
     .Set "UseRatioLimit", "0" 
     .Set "RatioLimit", "100" 
     .Set "MinStep", "0" 
     'MESHING METHOD 
     .SetMeshType "Unstr" 
     .Set "Method", "0" 
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "CurvatureOrder", "1" 
     .Set "CurvatureOrderPolicy", "off" 
     .Set "CurvRefinementControl", "NormalTolerance" 
     .Set "NormalTolerance", "22.5" 
     .Set "SrfMeshGradation", "2" 
     .Set "SrfMeshOptimization", "1" 
End With 
With MeshSettings 
     .SetMeshType "Unstr" 
     .Set "UseMaterials",  "1" 
End With 
With MeshSettings 
     .SetMeshType "Tet" 
     .Set "UseAnisoCurveRefinement", "1" 
     .Set "UseSameSrfAndVolMeshGradation", "1" 
     .Set "VolMeshGradation", "2" 
     .Set "VolMeshOptimization", "1" 
End With 
With MeshSettings 
     .SetMeshType "Unstr" 
     .Set "SmallFeatureSize", "0" 
     .Set "CoincidenceTolerance", "1e-006" 
     .Set "SelfIntersectionCheck", "1" 
     .Set "OptimizeForPlanarStructures", "0" 
End With 
With Mesh 
     .SetParallelMesherMode "Tet", "maximum" 
     .SetMaxParallelMesherThreads "Tet", "1" 
End With 

'@ set optimizer parameters

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Optimizer
  .SetMinMaxAuto "10" 
  .SetAlwaysStartFromCurrent "True" 
  .ResetParameterList
  .SelectParameter "r1", "False" 
  .SetParameterInit "1.0" 
  .SetParameterMin "0.9" 
  .SetParameterMax "1.1" 
  .SetParameterAnchors "5" 
End With 


'@ set optimizer settings

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Optimizer 
  .SetOptimizerType "Trust_Region" 
  .SetSimulationType "Frequency Domain Solver" 
  .SetAccuracy "0.01" 
  .SetNumRefinements "1" 
  .SetGenerationSize "32", "Genetic_Algorithm" 
  .SetGenerationSize "30", "Particle_Swarm" 
  .SetMaxIt "30", "Genetic_Algorithm" 
  .SetMaxIt "15", "Particle_Swarm" 
  .SetMaxEval "3000", "CMAES" 
  .SetUseMaxEval "True", "CMAES" 
  .SetSigma "0.2", "CMAES" 
  .SetSeed "1", "CMAES" 
  .SetSeed "1", "Genetic_Algorithm" 
  .SetSeed "1", "Particle_Swarm" 
  .SetSeed "1", "Nelder_Mead_Simplex" 
  .SetMaxEval "500", "Trust_Region" 
  .SetUseMaxEval "False", "Trust_Region" 
  .SetSigma "0.2", "Trust_Region" 
  .SetDomainAccuracy "0.01", "Trust_Region" 
  .SetFiniteDiff "0", "Trust_Region" 
  .SetMaxEval "250", "Nelder_Mead_Simplex" 
  .SetUseMaxEval "False", "Nelder_Mead_Simplex" 
  .SetUseInterpolation "No_Interpolation", "Genetic_Algorithm" 
  .SetUseInterpolation "No_Interpolation", "Particle_Swarm" 
  .SetInitialDistribution "Latin_Hyper_Cube", "Genetic_Algorithm" 
  .SetInitialDistribution "Latin_Hyper_Cube", "Particle_Swarm" 
  .SetInitialDistribution "Noisy_Latin_Hyper_Cube", "Nelder_Mead_Simplex" 
  .SetUsePreDefPointInInitDistribution "True", "Nelder_Mead_Simplex" 
  .SetUsePreDefPointInInitDistribution "True", "CMAES" 
  .SetGoalFunctionLevel "0.0", "Genetic_Algorithm" 
  .SetGoalFunctionLevel "0.0", "Particle_Swarm" 
  .SetGoalFunctionLevel "0.0", "Nelder_Mead_Simplex" 
  .SetMutaionRate "60", "Genetic_Algorithm" 
  .SetMinSimplexSize "1e-6" 
  .SetGoalSummaryType "Sum_All_Goals" 
  .SetUseDataOfPreviousCalculations "False" 
  .SetDataStorageStrategy "None" 
End With 


'@ define background

'[VERSION]2015.0|24.0.2|20150116[/VERSION]
With Background 
     .ResetBackground 
     .XminSpace "5" 
     .XmaxSpace "5" 
     .YminSpace "5" 
     .YmaxSpace "5" 
     .ZminSpace "5" 
     .ZmaxSpace "5" 
     .ApplyInAllDirections "True" 
End With 

With Material 
     .Reset 
     .FrqType "all"
     .Type "Normal"
     .MaterialUnit "Frequency", "Hz"
     .MaterialUnit "Geometry", "m"
     .MaterialUnit "Time", "s"
     .MaterialUnit "Temperature", "Kelvin"
     .Epsilon "1.0"
     .Mue "1.0"
     .Sigma "0.0"
     .TanD "0.0"
     .TanDFreq "0.0"
     .TanDGiven "False"
     .TanDModel "ConstSigma"
     .EnableUserConstTanDModelOrderEps "False"
     .ConstTanDModelOrderEps "1"
     .SetElParametricConductivity "False"
     .ReferenceCoordSystem "Global"
     .CoordSystemType "Cartesian"
     .SigmaM "0"
     .TanDM "0.0"
     .TanDMFreq "0.0"
     .TanDMGiven "False"
     .TanDMModel "ConstSigma"
     .EnableUserConstTanDModelOrderMue "False"
     .ConstTanDModelOrderMue "1"
     .SetMagParametricConductivity "False"
     .DispModelEps  "None"
     .DispModelMue "None"
     .DispersiveFittingSchemeEps "Nth Order"
     .MaximalOrderNthModelFitEps "10"
     .ErrorLimitNthModelFitEps "0.1"
     .UseOnlyDataInSimFreqRangeNthModelEps "False"
     .DispersiveFittingSchemeMue "Nth Order"
     .MaximalOrderNthModelFitMue "10"
     .ErrorLimitNthModelFitMue "0.1"
     .UseOnlyDataInSimFreqRangeNthModelMue "False"
     .UseGeneralDispersionEps "False"
     .UseGeneralDispersionMue "False"
     .NLAnisotropy "False"
     .NLAStackingFactor "1"
     .NLADirectionX "1"
     .NLADirectionY "0"
     .NLADirectionZ "0"
     .Rho "0.0"
     .ThermalType "Normal"
     .ThermalConductivity "0.0"
     .HeatCapacity "0.0"
     .MetabolicRate "0"
     .BloodFlow "0"
     .VoxelConvection "0"
     .MechanicsType "Unused"
     .Colour "0.6", "0.6", "0.6" 
     .Wireframe "False" 
     .Reflection "False" 
     .Allowoutline "True" 
     .Transparentoutline "False" 
     .Transparency "0" 
     .ChangeBackgroundMaterial
End With 


