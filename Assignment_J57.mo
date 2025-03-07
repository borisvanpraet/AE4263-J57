package Assignment_J57
  model Design
  inner BasicAeroEngines.Components.Environment designEnvironment annotation(
      Placement(transformation(origin = {70, 76}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CombustionChamberLHVSimple combustionChamberLHVSimple annotation(
      Placement(transformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Data.Compressors.SampleCompressor sampleCompressor annotation(
      Placement(transformation(origin = {-50, 32}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Data.Turbines.SampleTurbine sampleTurbine annotation(
      Placement(transformation(origin = {40, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.TimeTable fuelFlow annotation(
      Placement(transformation(origin = {-30, 70}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.AirIntake airIntake annotation(
      Placement(transformation(origin = {-96, -38}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CompressorMapsBetaLines LPC annotation(
      Placement(transformation(origin = {-62, -38}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CompressorBleed HPC annotation(
      Placement(transformation(origin = {-28, 0}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.TurbineMapsBetaLines HPT annotation(
      Placement(transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CooledTurbine cooledLPT annotation(
      Placement(transformation(origin = {54, -38}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.ShaftInertia HP_shaft annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.ShaftInertia LP_shaft annotation(
      Placement(transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.NozzleExhaust nozzleExhaust annotation(
      Placement(transformation(origin = {90, -38}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(fuelFlow.y, combustionChamberLHVSimple.fuelFlow) annotation(
      Line(points = {{-18, 70}, {0, 70}, {0, 40}}, color = {0, 0, 127}));
  connect(airIntake.outlet, LPC.inlet) annotation(
      Line(points = {{-90, -32}, {-68, -32}, {-68, -28}}, color = {170, 223, 255}));
  connect(LPC.outlet, HPC.inlet) annotation(
      Line(points = {{-56, -32}, {-34, -32}, {-34, 10}}, color = {170, 223, 255}));
  connect(HPC.outlet, combustionChamberLHVSimple.airInlet) annotation(
      Line(points = {{-22, 6}, {-10, 6}, {-10, 30}}, color = {170, 223, 255}));
  connect(combustionChamberLHVSimple.exhaust, HPT.inlet) annotation(
      Line(points = {{10, 30}, {20, 30}, {20, 6}}, color = {129, 170, 194}));
  connect(HPT.outlet, cooledLPT.inlet) annotation(
      Line(points = {{32, 10}, {48, 10}, {48, -32}}, color = {129, 170, 194}));
  connect(cooledLPT.outlet, nozzleExhaust.inlet) annotation(
      Line(points = {{60, -28}, {84, -28}, {84, -32}}, color = {129, 170, 194}));
  connect(HPT.shaft, HP_shaft.flange_b) annotation(
      Line(points = {{20, 0}, {10, 0}}));
  connect(HP_shaft.flange_a, HPC.shaft_b) annotation(
      Line(points = {{-10, 0}, {-22, 0}}));
  connect(cooledLPT.shaft, LP_shaft.flange_b) annotation(
      Line(points = {{48, -38}, {10, -38}, {10, -40}}));
  connect(LP_shaft.flange_a, LPC.shaft_b) annotation(
      Line(points = {{-10, -40}, {-56, -40}, {-56, -38}}));
    annotation(
      Icon(graphics = {Text(origin = {1, 1}, extent = {{-81, 79}, {81, -79}}, textString = "ADP")}));
end Design;
  annotation(
    Icon(graphics = {Text(origin = {3, 4}, extent = {{-63, 54}, {63, -54}}, textString = "J57")}),
  uses(BasicAeroEngines(version = "2.0.0"), Modelica(version = "4.0.0")));
end Assignment_J57;
