package HydroPowerControl "Simple models of hydropower plants"
  extends Modelica.Icons.Package;

  package Interfaces "Connectors"
    extends Modelica.Icons.InterfacesPackage;

    connector Flange "Connettore idraulico"
      Modelica.Units.SI.Pressure p "Pressure";
      flow Modelica.Units.SI.MassFlowRate w "Mass flow rate";
      annotation(
        Icon(graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid)}));
    end Flange;

    connector FlangeA "Hydraulic connector"
      extends Flange;
    end FlangeA;

    connector FlangeB "Hydraulic connector"
      extends Flange;
      annotation(
        Icon(graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid)}));
    end FlangeB;
    annotation(
      Diagram(graphics));
  end Interfaces;

  package Components "Component models"
    extends Modelica.Icons.Package;

    model FixedPressure "Ideal source at fixed pressure"
      parameter Modelica.Units.SI.Pressure p0 = 101325 "Fixed pressure";
      Interfaces.Flange flange annotation(
        Placement(transformation(extent = {{80, -10}, {100, 10}}, rotation = 0)));
    equation
      flange.p = p0;
      annotation(
        Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 255}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid)}));
    end FixedPressure;

    model FlowSource "Ideal prescribed flow source"
      Interfaces.Flange flange annotation(
        Placement(transformation(extent = {{80, -10}, {100, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput flowRate annotation(
        Placement(visible = true, transformation(origin = {-90, -2}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      flange.w = -flowRate;
      annotation(
        Icon(graphics = {Rectangle(origin = {9.56, -0.26}, lineColor = {0, 0, 255}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-70.12, 31.81}, {70.12, -31.81}}), Polygon(origin = {12, 2}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-30, 20}, {30, 0}, {-30, -20}, {0, 0}, {-30, 20}, {-30, 20}, {-30, 20}})}));
    end FlowSource;

    model Penstock "1D Penstock model with compressible fluid and elastic walls"
      parameter Modelica.Units.SI.Density rho = 999 "Fluid density";
      parameter Modelica.Units.SI.BulkModulus B = 21000e5 "Bulk Modulus of the fluid";
      parameter Modelica.Units.SI.Pressure E = 206e9 "Young modulus of the wall material";
      parameter Modelica.Units.SI.Length L "Penstock length";
      parameter Modelica.Units.SI.Length s "Wall thickness";
      parameter Modelica.Units.SI.Length D "Penstock diameter";
      parameter Integer N = 10 "Number of sections";
      parameter Modelica.Units.SI.Length hi = 0 "Inlet height";
      parameter Modelica.Units.SI.Length ho = 0 "Outlet height";
      parameter Boolean downstreamCapacitance = true "=true if capacitance is place at each volume outlet";
      parameter Real cf "Fanning friction factor";
      parameter Modelica.Units.SI.MassFlowRate wstart "Start value of mass flow rate";
      final parameter Modelica.Units.SI.Pressure K = E*s/D "Elastic coefficient of the pipe walls";
      final parameter Modelica.Units.SI.VelocityOfSound c = 1/sqrt(rho/B + rho/K) "Speed of sound";
      final parameter Modelica.Units.SI.Length z[N + 1] = linspace(hi, ho, N + 1) "Heights of individual volume boundaries";
      constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n "Acceleration of gravity";
      Modelica.Units.SI.Pressure p[N + 1] "Pressures at individual volumes boundaries";
      Modelica.Units.SI.Pressure pm[N] "Average pressures of volumes";
      Modelica.Units.SI.MassFlowRate w[N + 1] "Mass flow rates at the individual volume boundaries";
      Modelica.Units.SI.MassFlowRate wm[N] "Average volume flow rates";
      Modelica.Units.SI.Area A = Modelica.Constants.pi*D^2/4 "Cross section";
      Modelica.Units.SI.Length omega = Modelica.Constants.pi*D "Wet perimeter";
      Modelica.Units.SI.Length l = L/N "Length of individual volume";
      Interfaces.FlangeA inlet annotation(
        Placement(transformation(extent = {{-100, -10}, {-80, 10}}, rotation = 0)));
      Interfaces.FlangeB outlet annotation(
        Placement(transformation(extent = {{80, -10}, {100, 10}}, rotation = 0)));
    equation
// Mass balance for each volume
      for j in 1:N loop
        A*l/c^2*der(pm[j]) = w[j] - w[j + 1];
      end for;
// Momentum balance for each volume
      for i in 1:N loop
        l*der(wm[i]) + rho*A*g*(z[i + 1] - z[i]) + A*(p[i + 1] - p[i]) + cf*omega*l/(rho*A^2)*wm[i]*abs(wm[i]) = 0;
      end for;
// Component boundaries
      w[1] = inlet.w;
      w[N + 1] = -outlet.w;
      p[1] = inlet.p;
      p[N + 1] = outlet.p;
      for k in 1:N loop
        wm[k] = if downstreamCapacitance then w[k] else w[k + 1];
        pm[k] = if downstreamCapacitance then p[k + 1] else p[k];
      end for;
    initial equation
      for i in 1:N loop
        wm[i] = wstart;
        der(wm[i]) = 0;
      end for;
      annotation(
        Icon(graphics = {Rectangle(extent = {{-60, 20}, {60, -20}}, lineColor = {0, 0, 255}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid), Line(points = {{-80, 0}, {-60, 0}}, color = {0, 0, 255}), Line(points = {{60, 0}, {80, 0}}, color = {0, 0, 255})}),
        DymolaStoredErrors);
    end Penstock;

    model PeltonTurbine
      parameter Modelica.Units.SI.Pressure p_atm = 101325 "Atmospheric pressue";
      parameter Modelica.Units.SI.Density rho = 999 "Fluid density";
      parameter Modelica.Units.SI.Length r "Average rotor radius";
      parameter Modelica.Units.SI.PerUnit eta_n "Nominal efficiency";
      Modelica.Units.SI.MassFlowRate w "Mass flow rate";
      Modelica.Units.SI.Pressure dp "Inlet-outlet delta-p";
      Modelica.Units.SI.Velocity u "Nozzle outlet velocity";
      Modelica.Units.SI.Power P_hyd "Ideal hydraulic power";
      Modelica.Units.SI.Power P "Mechanical power";
      Modelica.Units.SI.PerUnit eta "Hydraulic efficiency";
      Modelica.Units.SI.PerUnit x;
      Modelica.Units.SI.AngularVelocity omega "Turbine angular velocity";
      Modelica.Units.SI.Torque tau "Torque applied to the mechanical shaft";
      Interfaces.FlangeA inlet "Inlet flange" annotation(
        Placement(transformation(extent = {{-100, 50}, {-80, 70}}, rotation = 0), iconTransformation(extent = {{-100, 40}, {-80, 60}})));
      Modelica.Blocks.Interfaces.RealInput An "Area of nozzle exhaust in m2" annotation(
        Placement(transformation(extent = {{-120, -14}, {-80, 28}}), iconTransformation(extent = {{-9.5, -9.5}, {9.5, 9.5}}, rotation = -90, origin = {-70.5, 90.5})));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft "Mechanical shaft" annotation(
        Placement(transformation(extent = {{-88, -50}, {-68, -30}}), iconTransformation(extent = {{-10, -10}, {10, 10}})));
    equation
      dp = inlet.p - p_atm;
      u = sqrt(2*dp/rho);
      w = rho*An*u;
      P_hyd = w*u^2/2;
      P = P_hyd*eta;
      omega*tau = P;
      eta = eta_n*(1 - (x - 1)^2);
      x = omega*r/u;
// Component boundaries
      w = inlet.w;
      omega = der(shaft.phi);
      tau = -shaft.tau;
      annotation(
        Icon(graphics = {Ellipse(extent = {{-52, 50}, {52, -50}}, lineColor = {0, 0, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-30, 30}, {30, -30}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-32, 62}, {-14, 48}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{-2, -50}, {12, -64}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{-52, 2}, {-66, -16}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{52, 8}, {66, 22}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{-40, 32}, {-62, 30}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{28, 42}, {24, 64}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{-36, -36}, {-32, -58}}, color = {0, 0, 255}, thickness = 1.25), Line(points = {{40, -32}, {64, -26}}, color = {0, 0, 255}, thickness = 1.25), Polygon(points = {{-70, 58}, {-70, 40}, {-50, 50}, {-70, 58}}, lineColor = {0, 0, 255}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid), Line(points = {{-80, 50}, {-72, 50}}, color = {0, 0, 255})}),
        DymolaStoredErrors);
    end PeltonTurbine;

    model GeneratorLoads "Model of ideal synchronous generator and loads"
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft "Generator shaft" annotation(
        Placement(visible = true, transformation(extent = {{-38, -8}, {-18, 12}}, rotation = 0), iconTransformation(extent = {{-110, -10}, {-90, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput loadReferencePower(unit = "W") "Reference value of active power" annotation(
        Placement(visible = true, transformation(extent = {{-42, -12}, {-2, 28}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      parameter Integer Np = 2 "Number of polar expansions of synchronous machine";
      parameter Modelica.Units.SI.PerUnit alpha = 1 "Exponent of frequency-power load curve";
      parameter Modelica.Units.SI.Frequency fn = 50 "Nominal grid frequency";
      final parameter Modelica.Units.SI.AngularVelocity omega_n = 2*pi*fn "Nominal electrical angular velocity";
      constant Modelica.Units.SI.PerUnit pi = Modelica.Constants.pi;
      Modelica.Units.SI.Power Pel "Active electrical power consumed by loads";
      Modelica.Units.SI.Torque tau "Turbine net torque";
      Modelica.Units.SI.AngularVelocity omega_m "Angular velocity of mechanical shaft";
      Modelica.Units.SI.AngularVelocity omega_e "Angular velocity of 3-phase AC system";
      Modelica.Blocks.Interfaces.RealOutput f(unit = "Hz") "Grid frequency in Hz" annotation(
        Placement(visible = true, transformation(extent = {{34, -10}, {54, 10}}, rotation = 0), iconTransformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
    equation
      Pel = omega_m*tau;
      omega_e = Np*omega_m;
      omega_e = f*2*pi;
      Pel = loadReferencePower*(omega_e/omega_n)^alpha;
// System boundary
      omega_m = der(shaft.phi);
      tau = shaft.tau;
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Rectangle(extent = {{-100, 102}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-100, 0}, {-70, 0}}, color = {0, 0, 0}), Ellipse(extent = {{-70, 20}, {-28, -20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{-30, 0}, {0, 0}}, color = {0, 0, 0}), Line(points = {{0, 40}, {0, -40}}, color = {0, 0, 0}), Line(points = {{0, 40}, {20, 40}}, color = {0, 0, 0}), Rectangle(extent = {{20, 48}, {44, 32}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{20, 8}, {44, -8}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{0, 0}, {20, 0}}, color = {0, 0, 0}), Rectangle(extent = {{20, -32}, {44, -48}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Line(points = {{0, -40}, {20, -40}}, color = {0, 0, 0})}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end GeneratorLoads;

    package Test
      extends Modelica.Icons.ExamplesPackage;

      model TestWavePropagation_N_10 "Test of wave propagation with N = 10"
        extends Modelica.Icons.Example;
        HydroPowerControl.Components.FlowSource flowSource annotation(
          Placement(visible = true, transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        HydroPowerControl.Components.Penstock pipe(D = 0.3, L = 1000, cf = 0, downstreamCapacitance = false, s = 0.15, wstart = 0) annotation(
          Placement(visible = true, transformation(origin = {18, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        HydroPowerControl.Components.FixedPressure sink annotation(
          Placement(visible = true, transformation(origin = {58, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Step step(height = 10, offset = 0, startTime = 0) annotation(
          Placement(visible = true, transformation(origin = {-58, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression delta_P_in_sim(y = pipe.inlet.p - sink.p0) annotation(
          Placement(visible = true, transformation(origin = {0, 80}, extent = {{-40, -10}, {40, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.RealExpression delta_P_exact(y = pipe.c/pipe.A*step.height*((if time > 0 then 1 else 0) + (if time > 2*pipe.L/pipe.c then -2 else 0))) annotation(
          Placement(visible = true, transformation(origin = {1, 48}, extent = {{-39, -10}, {39, 10}}, rotation = 0)));
      equation
        connect(step.y, flowSource.flowRate) annotation(
          Line(points = {{-47, 0}, {-33, 0}}, color = {0, 0, 127}));
        connect(flowSource.flange, pipe.inlet) annotation(
          Line(points = {{-15, 0}, {-1, 0}}, color = {0, 0, 255}));
        connect(pipe.outlet, sink.flange) annotation(
          Line(points = {{36, 0}, {50, 0}}, color = {0, 0, 255}));
        annotation(
          experiment(StartTime = -0.5, StopTime = 2.5, Tolerance = 1e-06, Interval = 0.002));
      end TestWavePropagation_N_10;

      model TestWavePropagation_N_40
        extends TestWavePropagation_N_10(pipe.N = 40);
        annotation(
          experiment(StartTime = -0.5, StopTime = 2.5, Tolerance = 1e-06, Interval = 0.002));
      end TestWavePropagation_N_40;
    end Test;
  end Components;

  package Systems "System models"
    extends Modelica.Icons.Package;

    model HydroPlant
      parameter Modelica.Units.SI.Power Pn (displayUnit = "MW")= 340e6 "Nominal active electrical power" annotation(Dialog(group="Generator data"));
      parameter Modelica.Units.SI.Frequency fn = 50 "Nominal grid frequency" annotation(Dialog(group="Generator data"));
      parameter Integer Np = 2 "Number of polar expansions of the generator" annotation(Dialog(group="Generator data"));
      parameter Modelica.Units.SI.Time Ta = 12 "Turbogenerator acceleration time constant" annotation(Dialog(group="Generator data"));
      parameter Modelica.Units.SI.Length h1 = 100 "Level of basin above penstock inlet" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.Length h2 = 900 "Penstock height difference" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.Length L = 1200 "Penstock length" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.Length D = 3 "Penstock internal diameter" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.Length s = 0.15 "Penstock wall thickness" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.PerUnit cf = 0.005 "Fanning friction factor" annotation(Dialog(group="Penstock data"));
      parameter Modelica.Units.SI.Length r = 0.892 "Pelton turbine radius" annotation(Dialog(group="Turbine data"));
      parameter Modelica.Units.SI.PerUnit eta_n = 0.88 "Pelton turbine nominal efficiency" annotation(Dialog(group="Turbine data"));
      parameter Modelica.Units.SI.MassFlowRate wstart = 40500 "Initial value of mass flow rate";
      parameter Modelica.Units.SI.Length rho = 999 "Density of cold water";
      parameter Modelica.Units.SI.Pressure p_atm = 101325 "Atmospheric pressure";
      
      constant Modelica.Units.SI.Acceleration g = 9.81 "Acceleration of gravity";
      final parameter Modelica.Units.SI.AngularVelocity omega_n = 2*pi*fn/Np "Nominal shaft angular velocity";
      final parameter Modelica.Units.SI.MomentOfInertia J = Pn*Ta/omega_n^2 "Turbogenerator moment of inertia";
      import Modelica.Constants.pi;
      Components.FixedPressure fixedPressure(p0 = p_atm + rho*g*h1) annotation(
        Placement(transformation(origin = {-78, 58}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Penstock penstock(L = L, s = s, D = D, cf = cf, hi = h2, ho = 0, wstart = wstart) annotation(
        Placement(transformation(origin = {-79, 9}, extent = {{-23, -23}, {23, 23}}, rotation = 270)));
      Components.PeltonTurbine peltonTurbine(r = r, eta_n = eta_n) annotation(
        Placement(transformation(extent = {{-38, -42}, {2, -2}})));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J = J, phi(fixed = true, start = 0), w(fixed = true, start = omega_n)) annotation(
        Placement(transformation(extent = {{12, -32}, {32, -12}})));
      Components.GeneratorLoads generatorLoads(Np = Np) annotation(
        Placement(transformation(extent = {{50, -32}, {70, -12}})));
      Modelica.Blocks.Interfaces.RealInput loadReferencePower "Reference value of active power" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {60, 28}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {0, 100})));
      Modelica.Blocks.Interfaces.RealInput An "Nozzle flow coefficient in m2" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {-32, 30}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = 0, origin = {-100, 0})));
      Modelica.Blocks.Interfaces.RealOutput f "Grid frequency" annotation(
        Placement(transformation(extent = {{82, -32}, {102, -12}}), iconTransformation(extent = {{82, -20}, {122, 20}})));
    equation
      connect(penstock.inlet, fixedPressure.flange) annotation(
        Line(points = {{-79, 29.7}, {-79, 49}, {-78, 49}}, color = {0, 0, 255}));
      connect(penstock.outlet, peltonTurbine.inlet) annotation(
        Line(points = {{-79, -11.7}, {-31.5, -11.7}, {-31.5, -12}, {-36, -12}}, color = {0, 0, 255}));
      connect(peltonTurbine.shaft, inertia.flange_a) annotation(
        Line(points = {{-18, -22}, {12, -22}}, color = {0, 0, 0}));
      connect(inertia.flange_b, generatorLoads.shaft) annotation(
        Line(points = {{32, -22}, {50, -22}}, color = {0, 0, 0}));
      connect(generatorLoads.loadReferencePower, loadReferencePower) annotation(
        Line(points = {{60, -12}, {60, 28}}, color = {0, 0, 127}));
      connect(peltonTurbine.An, An) annotation(
        Line(points = {{-32.1, -3.9}, {-32.1, 6.05}, {-32, 6.05}, {-32, 30}}, color = {0, 0, 127}));
      connect(generatorLoads.f, f) annotation(
        Line(points = {{70, -22}, {92, -22}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {102, -100}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-68, 64}, {78, -64}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, textString = "P")}));
    end HydroPlant;

    model TestStepResponseAn "Open loop response to step on An"
      extends Modelica.Icons.Example;
      HydroPlant hydroPlant annotation(
        Placement(transformation(extent = {{0, -20}, {40, 20}})));
      Modelica.Blocks.Sources.Step An(height = -0.029, offset = 0.29337, startTime = 30) annotation(
        Placement(transformation(extent = {{-68, -10}, {-48, 10}})));
      Modelica.Blocks.Sources.Step Pel(height = 0, offset = 340e6, startTime = 30) annotation(
        Placement(transformation(extent = {{-38, 54}, {-18, 74}})));
    equation
      connect(hydroPlant.An, An.y) annotation(
        Line(points = {{0, 0}, {-47, 0}}, color = {0, 0, 127}));
      connect(Pel.y, hydroPlant.loadReferencePower) annotation(
        Line(points = {{-17, 64}, {20, 64}, {20, 20}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)),
        experiment(StopTime = 60, StartTime = 0, Tolerance = 1e-06, Interval = 0.012));
    end TestStepResponseAn;

    model TestStepResponsePel "Open loop response to step on An"
      extends Modelica.Icons.Example;
      HydroPlant hydroPlant annotation(
        Placement(transformation(extent = {{0, -20}, {40, 20}})));
      Modelica.Blocks.Sources.Step An(height = 0, offset = 0.29337, startTime = 30) annotation(
        Placement(transformation(extent = {{-68, -10}, {-48, 10}})));
      Modelica.Blocks.Sources.Step Pel(height = -34e6, offset = 340e6, startTime = 30) annotation(
        Placement(transformation(extent = {{-38, 54}, {-18, 74}})));
    equation
      connect(hydroPlant.An, An.y) annotation(
        Line(points = {{0, 0}, {-47, 0}}, color = {0, 0, 127}));
      connect(Pel.y, hydroPlant.loadReferencePower) annotation(
        Line(points = {{-17, 64}, {20, 64}, {20, 20}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)),
        experiment(StopTime = 60, StartTime = 0, Tolerance = 1e-06, Interval = 0.12));
    end TestStepResponsePel;

    model ControlledSystemPI "Open loop response to step on An"
      extends Modelica.Icons.Example;
      HydroPowerControl.Systems.HydroPlant hydroPlant annotation(
        Placement(visible = true, transformation(extent = {{82, -20}, {122, 20}}, rotation = 0)));
      Modelica.Blocks.Continuous.PI PI(k = 0, T = 1, initType = Modelica.Blocks.Types.Init.InitialOutput) annotation(
        Placement(transformation(extent = {{-84, -10}, {-64, 10}})));
      Modelica.Blocks.Math.Feedback feedback annotation(
        Placement(transformation(extent = {{-116, -10}, {-96, 10}})));
      Modelica.Blocks.Sources.Constant fn(k = 50) annotation(
        Placement(transformation(extent = {{-146, -10}, {-126, 10}})));
      Modelica.Blocks.Math.Add add annotation(
        Placement(transformation(extent = {{-14, -4}, {6, 16}})));
      Modelica.Blocks.Sources.Constant An_(k = 0.2943) annotation(
        Placement(transformation(extent = {{-54, 16}, {-34, 36}})));
      Modelica.Blocks.Continuous.FirstOrder actuator(T = 0.2, initType = Modelica.Blocks.Types.Init.SteadyState) annotation(
        Placement(transformation(extent = {{38, -10}, {58, 10}})));
      Modelica.Blocks.Sources.TimeTable Pel(table = [0, 340e6; 30, 340e6; 30, 306e6; 200, 306e6]) annotation(
        Placement(visible = true, transformation(origin = {60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(PI.u, feedback.y) annotation(
        Line(points = {{-86, 0}, {-97, 0}}, color = {0, 0, 127}));
      connect(hydroPlant.f, feedback.u2) annotation(
        Line(points = {{122, 0}, {140, 0}, {140, -40}, {-106, -40}, {-106, -8}}, color = {0, 0, 127}));
      connect(fn.y, feedback.u1) annotation(
        Line(points = {{-125, 0}, {-114, 0}}, color = {0, 0, 127}));
      connect(PI.y, add.u2) annotation(
        Line(points = {{-63, 0}, {-16, 0}}, color = {0, 0, 127}));
      connect(An_.y, add.u1) annotation(
        Line(points = {{-33, 26}, {-28, 26}, {-28, 12}, {-16, 12}}, color = {0, 0, 127}));
      connect(actuator.y, hydroPlant.An) annotation(
        Line(points = {{59, 0}, {82, 0}}, color = {0, 0, 127}));
      connect(add.y, actuator.u) annotation(
        Line(points = {{7, 6}, {20, 6}, {20, 0}, {36, 0}}, color = {0, 0, 127}));
      connect(Pel.y, hydroPlant.loadReferencePower) annotation(
        Line(points = {{71, 50}, {102, 50}, {102, 20}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-160, -100}, {160, 100}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-160, -100}, {160, 100}})),
        experiment(StopTime = 60, StartTime = 0, Tolerance = 1e-06, Interval = 0.12));
    end ControlledSystemPI;
  end Systems;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false)),
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    uses(Modelica(version = "4.0.0")));
end HydroPowerControl;
