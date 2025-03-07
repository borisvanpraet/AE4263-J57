within BasicAeroEngines;
package BaseClasses "Base models"
  extends Modelica.Icons.BasesPackage;

  partial model BaseCompressor
    outer Components.Environment environment;
    package Air = Media.Air;
    replaceable parameter Types.TurbomachineData data "Compressor data";
    final parameter Modelica.Units.SI.Pressure P_E_nom(fixed=false)
      "Nominal/design inlet pressure";
    final parameter Modelica.Units.SI.Pressure P_L_nom(fixed=false)
      "Nominal/design outlet pressure";
    final parameter Modelica.Units.SI.Temperature T_E_nom(fixed=false)
      "Nominal inlet temperature";
    final parameter Modelica.Units.SI.PerUnit eta_nom=data.eta_nom
      "Nominal Isentropic efficiency";
    final parameter Modelica.Units.SI.MassFlowRate f_nom(fixed=false)
      "Nominal/design mass flow rate";
    final parameter Modelica.Units.SI.AngularVelocity omega_nom(fixed=false)
      "Nominal/design angular velocity";

    final parameter Modelica.Units.SI.SpecificHeatCapacity cp_nom=1000
      "Nominal cp";
    final parameter Modelica.Units.SI.PerUnit e=0.28 "(gamma-1)/gamma";
    final parameter Modelica.Units.SI.Torque tau_nom=cp_nom*data.T_E_nom*((data.P_L_nom
        /data.P_E_nom)^e - 1)/eta_nom*data.f_nom/data.omega_nom
      "Nominal torque";
    final parameter Modelica.Units.SI.Temperature T_L_start=data.T_E_nom*(data.P_L_nom
        /data.P_E_nom)^e "Start value of leaving temperature";
    final parameter Modelica.Units.SI.SpecificEnthalpy h_E_start=
        Air.specificEnthalpy_pTX(data.P_E_nom, data.T_E_nom)
      "Start value of entering enthalpy";
    final parameter Modelica.Units.SI.SpecificEnthalpy h_L_start=
        Air.specificEnthalpy_pTX(data.P_L_nom, T_L_start)
      "Start value of leaving enthalpy";

    Interfaces.AirPort inlet "Air inlet port" annotation (
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Interfaces.AirPort outlet "Air outlet port" annotation (
      Placement(visible = true, transformation(origin = {60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

    Modelica.Units.SI.MassFlowRate f_E "Entering mass flow rate";
    Modelica.Units.SI.MassFlowRate f_L "Leaving mass flow rate";
    Modelica.Units.SI.MassFlowRate f_B = 0 "Extra mass flow rate leaving before the outlet";
    Modelica.Units.SI.Pressure P_E(start=data.P_E_nom) "Entering pressure";
    Modelica.Units.SI.Pressure P_L(start=data.P_L_nom) "Leaving pressure";
    Modelica.Units.SI.SpecificEnthalpy h_E(start=h_E_start)
      "Entering specific enthalpy";
    Modelica.Units.SI.SpecificEnthalpy h_L(start=h_L_start)
      "Leaving specific enthalpy";
    Modelica.Units.SI.SpecificEnthalpy h_iso(start=h_L_start)
      "Isentropic enthalpy at outlet pressure";
    Modelica.Units.SI.Power W_B = 0 "Extra power flow leaving before the outlet";
    Modelica.Units.SI.Temperature T_E(start=data.T_E_nom) "Inlet temperature";
    Modelica.Units.SI.Temperature T_L(start=data.T_E_nom*(data.P_L_nom/data.P_E_nom)
          ^e) "Outlet temperature";
    Modelica.Units.SI.PerUnit PR "Pressure ratio";
    Modelica.Units.SI.PerUnit eta "Isentropic efficiency";
    Modelica.Units.SI.Power W "Mechanical power input to the compressor";
    Modelica.Units.SI.AngularVelocity omega "Angular velocity of shaft";
    Modelica.Units.SI.Torque tau(start=tau_nom)
      "Torque applied onto the shaft (positive)";
    Modelica.Units.SI.PerUnit PR_n "Normalized pressure ratio";
    Modelica.Units.SI.PerUnit Phi_n "Normalized flow number";
    //Real Phi_c_DP "Flow number @ design point";
    Real Phi_c "Flow number";
    Modelica.Units.SI.PerUnit N_n(start=1) "Normalized corrected speed";
    Modelica.Units.SI.PerUnit eta_n "Normalized Isentropic efficiency";
    Air.BaseProperties propIn "Fluid properties at the inlet";
    Air.BaseProperties propOut "Fluid properties at the outlet";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft_a annotation (
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Interfaces.Flange_b shaft_b annotation (
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  equation
// Balance equations
    f_E = f_L + f_B "Mass balance";
    W = f_L*h_L - f_E*h_E + W_B "Energy balance";
// Definition of auxiliary quantities
    W = omega * tau;
    PR = P_L / P_E;
    PR_n = PR / (P_L_nom / P_E_nom);
    eta = (h_iso - h_E) / (h_L - h_E);
    eta_n = eta / eta_nom;
    Phi_c = f_E * sqrt(T_E) / P_E;
//Phi_c_DP = (f_nom*sqrt(T_E_nom)/P_E_nom);
    Phi_n = Phi_c / (f_nom * sqrt(T_E_nom) / P_E_nom);
    N_n = omega / sqrt(T_E) / (omega_nom / sqrt(T_E_nom));
// Fluid properties
    propIn.p = P_E;
    propIn.h = h_E;
    propOut.p = P_L;
    propOut.h = h_L;
    T_E = propIn.T;
    T_L = propOut.T;
    h_iso = Air.isentropicEnthalpy(P_L, propIn.state);
// Boundary conditions
    f_E = inlet.f;
    f_L = -outlet.f;
    P_E = inlet.P;
    P_L = outlet.P;
    h_E = inStream(inlet.h_L);
    h_L = outlet.h_L;
    inlet.h_L = 0 "Not used, no flow reversal";
    shaft_a.phi = shaft_b.phi;
    omega = der(shaft_a.phi);
    tau = shaft_a.tau + shaft_b.tau;
    annotation (
      Icon(graphics={  Polygon(fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, points = {{-60, 100}, {-60, -100}, {60, -60}, {60, 60}, {60, 60}, {-60, 100}})}));
  end BaseCompressor;

  partial model BaseTurbine
    package ExhaustGas = Media.ExhaustGas;
    outer Components.Environment environment;

    constant Modelica.Units.SI.PerUnit gamma=1.4
      "Approximated cp/cv ratio for start value computation";

    final parameter Modelica.Units.SI.Pressure P_E_nom(fixed=false)
      "Nominal/design inlet pressure";
    final parameter Modelica.Units.SI.Pressure P_L_nom(fixed=false)
      "Nominal/design outlet pressure";
    final parameter Modelica.Units.SI.Temperature T_E_nom(fixed=false)
      "Nominal inlet temperature";
    final parameter Modelica.Units.SI.PerUnit eta_nom=data.eta_nom
      "Nominal Isentropic efficiency";
    final parameter Modelica.Units.SI.MassFlowRate f_nom(fixed=false)
      "Nominal/design mass flow rate";
    final parameter Modelica.Units.SI.AngularVelocity omega_nom(fixed=false)
      "Nominal/design angular velocity";

    replaceable parameter Types.TurbomachineData data "Turbine data";
    final parameter Modelica.Units.SI.SpecificHeatCapacity cp_nom=1000
      "Nominal cp";
    final parameter Modelica.Units.SI.PerUnit e=0.28 "(gamma-1)/gamma";
    final parameter Modelica.Units.SI.Torque tau_nom=cp_nom*data.T_E_nom*(1 - (
        data.P_L_nom/data.P_E_nom)^e)*eta_nom*data.f_nom/data.omega_nom
      "Nominal torque";
    final parameter Modelica.Units.SI.Temperature T_L_start=data.T_E_nom*(data.P_L_nom
        /data.P_E_nom)^e "Start value of leaving temperature";
    final parameter Modelica.Units.SI.SpecificEnthalpy h_E_start=
        ExhaustGas.specificEnthalpy_pTX(data.P_E_nom, data.T_E_nom)
      "Start value of entering enthalpy";
    final parameter Modelica.Units.SI.SpecificEnthalpy h_L_start=
        ExhaustGas.specificEnthalpy_pTX(data.P_L_nom, T_L_start)
      "Start value of leaving enthalpy";
    Interfaces.ExhaustPort inlet annotation (
      Placement(visible = true, transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Interfaces.ExhaustPort outlet annotation (
      Placement(visible = true, transformation(origin = {60, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Units.SI.MassFlowRate f_E "Entering mass flow rate";
    Modelica.Units.SI.MassFlowRate f_L "Leaving mass flow rate";
    Modelica.Units.SI.MassFlowRate f_B = 0 "Extra mass flow rate entering before the outlet";
    Modelica.Units.SI.Pressure P_E(start=data.P_E_nom) "Entering pressure";
    Modelica.Units.SI.Pressure P_L(start=data.P_L_nom) "Leaving pressure";
    Modelica.Units.SI.MassFraction X_E[ExhaustGas.nX]
      "Entering gas composition";
    Modelica.Units.SI.MassFraction X_L[ExhaustGas.nX] = X_E "Leaving gas composition; X_L = X_E Single components mass balances";
    Modelica.Units.SI.SpecificEnthalpy h_E(start=h_E_start)
      "Entering specific enthalpy";
    Modelica.Units.SI.SpecificEnthalpy h_L(start=h_L_start)
      "Leaving specific enthaly";
    Modelica.Units.SI.SpecificEnthalpy h_out = h_L
      "Specific enthaly for the calculation of isentropic efficiency";
    Modelica.Units.SI.SpecificEnthalpy h_iso
      "Isentropic enthalpy at outlet pressure";
    Modelica.Units.SI.Power W_B = 0 "Extra power flow leaving before the outlet";
    Modelica.Units.SI.Temperature T_E(start=data.T_E_nom) "Inlet temperature";
    Modelica.Units.SI.Temperature T_L(start=T_L_start)
                                                      "Outlet temperature";
    Modelica.Units.SI.PerUnit PR(start=data.P_E_nom/data.P_L_nom)
      "Pressure ratio";
    Modelica.Units.SI.PerUnit eta(start=0.8) "Isentropic efficiency";
    Modelica.Units.SI.Power W "Mechanical power";
    Modelica.Units.SI.AngularVelocity omega "Angular velocity of shaft";
    Modelica.Units.SI.Torque tau(start=tau_nom)
      "Torque applied onto the shaft (positive)";
    Modelica.Units.SI.PerUnit PR_n(start=1) "Normalized pressure ratio";
    Modelica.Units.SI.PerUnit Phi_n(start=1) "Normalized flow number";
    Modelica.Units.SI.PerUnit N_n(start=1) "Normalized corrected speed";
    Modelica.Units.SI.PerUnit eta_n "Normalized Isentropic efficiency";
    ExhaustGas.BaseProperties propIn "Fluid properties at the inlet";
    ExhaustGas.BaseProperties propOut "Fluid properties at the outlet";
    parameter Real eta_mech = 0.99 "Mechanical efficiency";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft annotation (
      Placement(visible = true, transformation(origin = {58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
// Balance equations
    f_E + f_B = f_L "Total mass balance";
    W = f_E*h_E - f_L*h_L + W_B "Energy balance";
// Definition of auxiliary quantities
    W = omega*tau;
    PR = P_E/P_L;
    PR_n = PR / (P_E_nom / P_L_nom);
    eta = (h_E - h_out) / (h_E - h_iso);
    eta_n = eta / eta_nom;
    Phi_n = f_E * sqrt(T_E) / P_E / (f_nom * sqrt(T_E_nom) / P_E_nom);
    N_n = omega/sqrt(T_E) / (omega_nom/sqrt(T_E_nom));
// Fluid properties
    propIn.p = P_E;
    propIn.h = h_E;
    propIn.X = X_E;
    propOut.p = P_L;
    propOut.h = h_L;
    propOut.X = X_E;
    T_E = propIn.T;
    T_L = propOut.T;
    h_iso = ExhaustGas.isentropicEnthalpy(P_L, propIn.state);
// Boundary conditions
    f_E = inlet.f;
    f_L = -outlet.f;
    P_E = inlet.P;
    P_L = outlet.P;
    h_E = inStream(inlet.h_L);
    h_L = outlet.h_L;
    X_E = inStream(inlet.X_L);
    X_L = outlet.X_L;
    inlet.h_L = 0 "Not used, no flow reversal";
    inlet.X_L = ExhaustGas.reference_X;
    omega = der(shaft.phi);
    tau = -shaft.tau / eta_mech;
    annotation (
      Icon(graphics={  Polygon(fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, points = {{-60, 60}, {-60, -60}, {60, -100}, {60, 100}, {60, 100}, {-60, 60}})}, coordinateSystem(initialScale = 0.1)));
  end BaseTurbine;
end BaseClasses;
