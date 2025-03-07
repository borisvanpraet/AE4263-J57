package Assignment_J57
  model Design
  inner BasicAeroEngines.Components.Environment designEnvironment(onDesignInit = true, altitude = 0, airspeed = 0) annotation(
      Placement(transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}})));
    //Intake
    BasicAeroEngines.Components.AirIntake airIntake annotation(
      Placement(transformation(origin = {-90, -38}, extent = {{-16, -16}, {16, 16}})));
    
    //Compressors
    BasicAeroEngines.Components.CompressorMapsBetaLines LPC(data = LPC_map) annotation(
    Placement(transformation(origin = {-60, -38}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.CompressorMapsBetaLines HPC(data = HPC_map) annotation(
      Placement(transformation(origin = {-32, -14}, extent = {{-10, -10}, {10, 10}})));
  
    //Combustor
    BasicAeroEngines.Components.CombustionChamberLHV combustionChamberLHV(LHV(displayUnit = "MJ/kg") = 4.296e7, V = 0.1, P_start(displayUnit = "kPa") = 1.120528e6, T_start(displayUnit = "K") = 598.4, steadyStateInit = true, ZC = 1, ZH = 1.9167, eta_comb = 0.9945, P_Loss = 5.5) annotation(
      Placement(transformation(origin = {0, 14}, extent = {{-10, -10}, {10, 10}})));
    
    //Turbines
    BasicAeroEngines.Components.CooledTurbine cooledHPT(Nstages = 1, Nbleed = 1, Nstages_Bleeds = {1}, Xi_cool_stat = {0.05}, Xi_cool_rot = {0.05}, eta_mech = 0.982, data = HPT_map) annotation(
        Placement(transformation(origin = {32, -14}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.TurbineMapsBetaLines LPT(eta_mech = 1, data = LPT_map) annotation(
      Placement(transformation(origin = {60, -38}, extent = {{-10, -10}, {10, 10}})));
    
    //Exhaust
    BasicAeroEngines.Components.NozzleExhaust nozzleExhaust(f_nom = 54.28, A_fixed = 0.2375) annotation(
    Placement(transformation(origin = {90, -38}, extent = {{-16, -16}, {16, 16}})));
  
    //Shafts
    BasicAeroEngines.Components.ShaftInertia HP_shaft(J = 74, omega_nom(displayUnit = "rpm") = 967.0869385300581) annotation(
      Placement(transformation(origin = {0, -14}, extent = {{-10, -10}, {10, 10}})));
    BasicAeroEngines.Components.ShaftInertia LP_shaft(J = 380, omega_nom(displayUnit = "rpm") = 623.606141737574) annotation(
      Placement(transformation(origin = {0, -38}, extent = {{-10, -10}, {10, 10}})));
    
    //Maps
    parameter BasicAeroEngines.Data.Compressors.GSP_LPC LPC_map(P_L_nom(displayUnit = "kPa") = 326684, T_E_nom(displayUnit = "K"), P_E_nom(displayUnit = "kPa"), omega_nom(displayUnit = "rpm"), eta_nom = 0.82) annotation(
      Placement(transformation(origin = {-58, -60}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Compressors.GSP_HPC HPC_map(P_L_nom(displayUnit = "kPa") = 1.120528e6, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.85) annotation(
      Placement(transformation(origin = {-40, 12}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Turbines.GSP_HPT_BetaLines HPT_map(P_L_nom(displayUnit = "kPa") = 462389, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.865) annotation(
      Placement(transformation(origin = {44, 12}, extent = {{-10, -10}, {10, 10}})));
    parameter BasicAeroEngines.Data.Turbines.GSP_LPT_BetaLines LPT_map(P_L_nom(displayUnit = "kPa") = 237732, P_E_nom(displayUnit = "kPa"), T_E_nom(displayUnit = "K"), omega_nom(displayUnit = "rpm"), eta_nom = 0.882) annotation(
      Placement(transformation(origin = {60, -60}, extent = {{-10, -10}, {10, 10}})));
  
    //Fuel input
    Modelica.Blocks.Sources.Constant fuelFlow(k = 1.03684) annotation(
      Placement(transformation(origin = {-18, 42}, extent = {{-10, -10}, {10, 10}})));
    
    //Expressions
    Modelica.Blocks.Sources.RealExpression netThrust(y = nozzleExhaust.thrust - airIntake.drag) annotation(
      Placement(transformation(origin = {142, -18}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression thermalEfficiency(y = (nozzleExhaust.W - airIntake.W)/combustionChamberLHV.Q) annotation(
      Placement(transformation(origin = {142, -38}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression overalEfficiency(y = (netThrust.y*designEnvironment.v)/combustionChamberLHV.Q) annotation(
      Placement(transformation(origin = {142, -60}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(LPT.outlet, nozzleExhaust.inlet) annotation(
      Line(points = {{66, -28}, {80, -28}}, color = {129, 170, 194}));
    connect(LPT.shaft, LP_shaft.flange_b) annotation(
      Line(points = {{54, -38}, {10, -38}}));
    connect(HP_shaft.flange_b, cooledHPT.shaft) annotation(
      Line(points = {{10, -14}, {26, -14}}));
    connect(combustionChamberLHV.exhaust, cooledHPT.inlet) annotation(
      Line(points = {{10, 14}, {26, 14}, {26, -8}}, color = {129, 170, 194}));
    connect(fuelFlow.y, combustionChamberLHV.fuelFlow) annotation(
      Line(points = {{-6, 42}, {0, 42}, {0, 24}}, color = {0, 0, 127}));
    connect(cooledHPT.outlet, LPT.inlet) annotation(
      Line(points = {{38, -4}, {54, -4}, {54, -32}}, color = {129, 170, 194}));
    connect(airIntake.outlet, LPC.inlet) annotation(
      Line(points = {{-80, -28}, {-66, -28}}, color = {170, 223, 255}));
  connect(HPC.shaft_b, HP_shaft.flange_a) annotation(
      Line(points = {{-26, -14}, {-10, -14}}));
  connect(HPC.outlet, combustionChamberLHV.airInlet) annotation(
      Line(points = {{-26, -8}, {-26, 14}, {-10, 14}}, color = {170, 223, 255}));
  connect(HPC.inlet, LPC.outlet) annotation(
      Line(points = {{-38, -4}, {-54, -4}, {-54, -32}}, color = {170, 223, 255}));
  connect(cooledHPT.Bl_port[1], HPC.outlet) annotation(
      Line(points = {{32, -6}, {32, -28}, {-26, -28}, {-26, -8}}, color = {170, 223, 255}));
  connect(LP_shaft.flange_a, LPC.shaft_a) annotation(
      Line(points = {{-10, -38}, {-66, -38}}));
    //annotation(
    //  __OpenModelica_commandLineOptions = "--tearingStrictness=veryStrict");
  annotation(
      Icon(graphics = {Text(origin = {1, 1}, extent = {{-81, 79}, {81, -79}}, textString = "ADP")}),
      Diagram,
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02));
  end Design;
  annotation(
    Icon(graphics = {Text(origin = {3, 4}, extent = {{-63, 54}, {63, -54}}, textString = "J57")}),
    uses(BasicAeroEngines(version = "2.0.0"), Modelica(version = "4.0.0")));
end Assignment_J57;
