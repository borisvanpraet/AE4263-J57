package Assignment_J57
  model Design
  inner BasicAeroEngines.Components.Environment designEnvironment annotation(
      Placement(transformation(origin = {70, 76}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.AirIntake airIntake annotation(
      Placement(transformation(origin = {-94, 0}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CompressorBleed compressorBleed annotation(
      Placement(transformation(origin = {-52, -12}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CombustionChamberLHVSimple combustionChamberLHVSimple annotation(
      Placement(transformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.CooledTurbine cooledTurbine annotation(
      Placement(transformation(origin = {44, -12}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.NozzleExhaust nozzleExhaust annotation(
      Placement(transformation(origin = {94, 0}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Components.ShaftInertia shaftInertia annotation(
      Placement(transformation(origin = {-2, -22}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Data.Compressors.SampleCompressor sampleCompressor annotation(
      Placement(transformation(origin = {-50, 32}, extent = {{-10, -10}, {10, 10}})));
  BasicAeroEngines.Data.Turbines.SampleTurbine sampleTurbine annotation(
      Placement(transformation(origin = {40, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.TimeTable fuelFlow annotation(
      Placement(transformation(origin = {-30, 70}, extent = {{-10, -10}, {10, 10}})));
  equation
  connect(cooledTurbine.shaft, shaftInertia.flange_b) annotation(
      Line(points = {{38, -12}, {23, -12}, {23, -22}, {8, -22}}));
  connect(shaftInertia.flange_a, compressorBleed.shaft_b) annotation(
      Line(points = {{-12, -22}, {-29, -22}, {-29, -12}, {-46, -12}}));
  connect(compressorBleed.outlet, combustionChamberLHVSimple.airInlet) annotation(
      Line(points = {{-46, -6}, {-29, -6}, {-29, 12}, {-10, 12}, {-10, 30}}, color = {170, 223, 255}));
  connect(combustionChamberLHVSimple.exhaust, cooledTurbine.inlet) annotation(
      Line(points = {{10, 30}, {23, 30}, {23, 12}, {38, 12}, {38, -6}}, color = {129, 170, 194}));
  connect(cooledTurbine.outlet, nozzleExhaust.inlet) annotation(
      Line(points = {{50, -2}, {88, -2}, {88, 6}}, color = {129, 170, 194}));
  connect(airIntake.outlet, compressorBleed.inlet) annotation(
      Line(points = {{-88, 6}, {-58, 6}, {-58, -2}}, color = {170, 223, 255}));
  connect(fuelFlow.y, combustionChamberLHVSimple.fuelFlow) annotation(
      Line(points = {{-18, 70}, {0, 70}, {0, 40}}, color = {0, 0, 127}));
  annotation(
      Icon(graphics = {Text(origin = {1, 1}, extent = {{-81, 79}, {81, -79}}, textString = "ADP")}));
end Design;
  annotation(
    Icon(graphics = {Text(origin = {3, 4}, extent = {{-63, 54}, {63, -54}}, textString = "J57")}),
  uses(BasicAeroEngines(version = "2.0.0"), Modelica(version = "4.0.0")));
end Assignment_J57;
