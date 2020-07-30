package projectWaterHeater
  /*
    Authors: Kristian Dyb Strand
    */

  model SimSimpleWaterHeater
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Simulation of model "ModSimpleWaterHeater"
      */
    // Instatiating model
    ModSimpleWaterHeater swh;
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real up = 1 "Output to heater, 0.0-1.0";
    parameter Real uv = 1;
    // Input variables
    //Real uv "Output to shunt valve, 0.0-1.0";
    // Output variable
    Real Te(unit = "°C") "Effluent temperature";
  equation
    swh.Ti = Ti;
    swh.Ta = Ta;
    swh.Vd = Vd;
//uv = if time < 3600 then 0 else 1;
    swh.uv = uv;
    swh.up = up;
    Te = swh.T;
  end SimSimpleWaterHeater;

  model SimSlicedWaterHeater1
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 28.10.19
      Purpose: Simulation of model "ModSlicedWaterHeater"
      */
    // Instatiating model
    //ModSlicedWaterHeater swh;
    ModSlicedWaterHeater2 swh2;
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real up = 1 "Output to heater, 0.0-1.0";
    parameter Real uv = 1 "Output to shunt valve, 0.0-1.0";
    // Input variables
    // Output variable
    Real Te(unit = "°C") "Effluent temperature";
    Real T[swh2.N] "Tank temperature profile";
  equation
/*
      swh.Ti = Ti;
      swh.Ta = Ta;
      swh.Vd = Vd;
      swh.uv = uv;
      swh.up = up;
    */
    swh2.Ti = Ti;
    swh2.Ta = Ta;
    swh2.Vd = Vd;
    swh2.uv = uv;
    swh2.up = up;
    Te = swh2.T[swh2.N];
    T[:] = swh2.T[:];
  end SimSlicedWaterHeater1;

  model SimSlicedWaterHeater2
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Simulation of model "ModSlicedWaterHeater"
      */
    // Instatiating model
    ModSlicedWaterHeater swh;
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real up = 1 "Output to heater, 0.0-1.0";
    // Input variables
    Real uv "Output to shunt valve, 0.0-1.0";
    // Output variable
    Real Te(unit = "°C") "Effluent temperature";
    Real T[swh.N] "Tank temperature profile";
  equation
    swh.Ti = Ti;
    swh.Ta = Ta;
    swh.Vd = Vd;
    uv = if time < 3600 then 1 else (7200 - (time - 3600)) / 7200;
    swh.uv = uv;
    swh.up = up;
    Te = swh.T[swh.N];
    T[:] = swh.T[:];
  end SimSlicedWaterHeater2;
  
  model SimSimpleWaterHeaterMatrix
  /*
    Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
    Revision: 1
    Date: 25.10.19
    Purpose: Simulation of model "ModSimpleWaterHeater"
    */
    // Instatiating model
    ModSimpleWaterHeaterMatrix swh;
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real up = 1 "Output to heater, 0.0-1.0";
    parameter Real uv = 1;
    // Input variables
    // Output variable
    Real Te(unit = "°C") "Effluent temperature";
  equation
    swh.Ti = Ti;
    swh.Ta = Ta;
    swh.Vd = Vd;
//uv = if time < 3600 then 0 else 1;
    swh.uv = uv;
    swh.up = up;
    Te = swh.T;
  end SimSimpleWaterHeaterMatrix;


  model SimShuntValve
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Simulation of model "ModShuntValve"
      */
    // Instantiate model
    ModShuntValve sv;
    // Parameters
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real Ti2(unit = "°C") = 27 "Influent temperature from loop";
    parameter Real Te_sp(unit = "°C") = 32 "Influent temperature from loop";
    // Input variables
    Real uv "Output to shunt valve, 0.0-1.0";
    Real Ti1(unit = "°C") "Influent temperature from tank";
    // Output variable
    Real Te "Effluent temperature";
  equation
    Te_sp = uv * Ti1 + (1 - uv) * Ti2;
    sv.uv = min(max(uv, 0), 1);
    sv.Vd = Vd;
    Ti1 = if time < 1 then 25 else 40;
    sv.Ti1 = Ti1;
    sv.Ti2 = Ti2;
    Te = sv.Te;
  end SimShuntValve;

  model SimTankAndValve
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Simulation of models "ModSlicedWaterHeater and "ModShuntValve" connected together in a system
      */
    // Instatiating model
    ModShuntValve sv;
    ModSlicedWaterHeater swh;
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real Te_sp(unit = "°C") = 32 "Influent temperature from loop";
    // Input variables
    Real uv "Output to shunt valve, 0.0-1.0";
    Real up "Output to heater, 0.0-1.0";
    // Output variables
    Real Te "Effluent temperature from system)";
    Real Tt "Effluent temperature from tank";
  equation
    up = 1 - (swh.T[swh.HeaterPosition] - 40) / 3;
    Tt = swh.T[swh.N];
    Te_sp = uv * Tt + (1 - uv) * Ti;
    sv.uv = min(max(uv, 0), 1);
    swh.uv = min(max(uv, 0), 1);
    swh.up = min(max(up, 0), 1);
    sv.Vd = Vd;
    swh.Vd = Vd;
    sv.Ti1 = Tt;
    sv.Ti2 = Ti;
    swh.Ti = Ti;
    swh.Ta = Ta;
    Te = sv.Te;
  end SimTankAndValve;

  model ModShuntValve
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Model for a Shuntvalve with two inlet streams and one outlet
      */
    // Constants
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    // Parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    // Auxiliary variables
    Real mdi1(unit = "kg/s") "Influent flow from loop";
    Real mdi2(unit = "kg/s") "Influent flow from tank";
    Real mde(unit = "kg/s") "Effluent flow";
    Real hdi1(unit = "W") "Influent enthalpy flow from tank";
    Real hdi2(unit = "W") "Influent enthalpy flow from loop";
    Real hde(unit = "W") "Efluent enthalpy flow";
    Real md(unit = "kg/s") "Total mass flow";
    // Input parameters
    input Real uv "Valve position";
    input Real Vd(unit = "l/min") "Total volum flow";
    input Real Ti1(unit = "°C") "Influent temperature from tank";
    input Real Ti2(unit = "°C") "Influent temperature from loop";
    // Output variable
    output Real Te(unit = "°C") "Effluent temperature";
  equation
// Algebraic equations
    mde = mdi1 + mdi2 "Mass balance";
    mdi1 = md * uv;
    mdi2 = md * (1 - uv);
    hde = hdi1 + hdi2 "Enthalpi balance";
    hde = mde * cp * (Te - Tref);
    hdi1 = mdi1 * cp * (Ti1 - Tref);
    hdi2 = mdi2 * cp * (Ti2 - Tref);
    md = rho * (Vd / (l2m3 * m2s));
  end ModShuntValve;

  model ModSimpleWaterHeater
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 25.10.19
      Purpose: Model of a perfectly mixed water heater
      */
    // Constants
    constant Real PI = 3.14;
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    // Parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    parameter Real d(unit = "m") = 0.5 "Tank diameter";
    parameter Real h(unit = "m") = 1.5 "Tank height";
    parameter Real A(unit = "m²") = 0.25 * PI * d ^ 2 "Tank bottom and top area";
    parameter Real V(unit = "m³") = A * h "Tank volume";
    parameter Real As(unit = "m²") = 2 * A + PI * d * h "Tank surface area";
    parameter Real m(unit = "kg") = rho * V "Water mass in tank";
    parameter Real ha(unit = "W/m²K") = 3 "Heat transfer, Air";
    parameter Real hw_free(unit = "W/m²K") = 50 "Heat transfer, Water, Free convection";
    parameter Real hw_forced(unit = "W/m²K") = 1000 "Heat transfer, Water, Forced convection";
    parameter Real d1k1(unit = "m²K/W") = 1 / 0.5 "Typical value for 5cm insulator";
    parameter Real P0(unit = "W") = 15e3 "Maximum power, Heating element";
    // Initial Conditions
    parameter Real T0(unit = "°C") = 25 "Initial temperature";
    parameter Real H0(unit = "J/kg") = cp * (T0 - Tref) "Intial specific enthalpy";
    parameter Real U0(unit = "J") = m * H0 "Initial internal energy";
    // State variables
    Real U(start = U0, fixed = true, unit = "J") "Internal energy";
    // Auxiliary Variables
    Real T(unit = "°C") "Temperature in tank";
    Real Hm(unit = "J/kg") "Specific enthalpy in tank";
    Real Hdi(unit = "W") "Influent enthalpy flow";
    Real Hde(unit = "W") "Efluent enthalpy flow";
    Real Hmi(unit = "J/kg") "Specific enthalpy in influent flow";
    Real Qd(unit = "W") "Total heat entering the system";
    Real Qdel(unit = "W") "Heat energy from electric heater";
    Real Qda(unit = "W") "Heat energy from surroundings";
    Real md(unit = "kg/s") "Mass flow through tank";
    Real Uloss(unit = "W/m²K") "Thermal heat transfer unit to surroundings";
    Real hw(unit = "W/m²K") "Heat transfer, Water";
    // Inputs
    input Real Ti(unit = "°C") "Influent temperatur";
    input Real Ta(unit = "°C") "Ambient temperature";
    input Real Vd(unit = "l/min") "Total volum flow";
    input Real up "Output to heater, 0.0-1.0";
    input Real uv "Output to shunt valve, 0.0-1.0";
    // Outputs
    //output Real T;
  equation
// Differential equations
    der(U) = Hdi - Hde + Qd;
// Algebraic equations
    U = m * Hm;
    Hdi = md * Hmi;
    Hde = md * Hm;
    Hmi = cp * (Ti - Tref);
    Hm = cp * (T - Tref);
    Qd = Qdel + Qda;
    Qdel = P0 * up;
    Qda = Uloss * As * (Ta - T);
    Uloss = 1 / (1 / ha + d1k1 + 1 / hw);
    hw = if uv < 1e-7 then hw_free else hw_forced;
    md = rho * (Vd / (l2m3 * m2s)) * uv;
  end ModSimpleWaterHeater;
  
  model ModSimpleWaterHeaterMatrix
  /*
    Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
    Revision: 1
    Date: 25.10.19
    Purpose: Model of a perfectly mixed water heater
    */
    // Constants
    constant Real PI = 3.14;
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    // Parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    parameter Real d(unit = "m") = 0.5 "Tank diameter";
    parameter Real h(unit = "m") = 1.5 "Tank height";
    parameter Real A(unit = "m²") = 0.25 * PI * d ^ 2 "Tank bottom and top area";
    parameter Real V(unit = "m³") = A * h "Tank volume";
    parameter Real As(unit = "m²") = 2 * A + PI * d * h "Tank surface area";
    parameter Real m(unit = "kg") = rho * V "Water mass in tank";
    parameter Real ha(unit = "W/m²K") = 3 "Heat transfer, Air";
    parameter Real hw_free(unit = "W/m²K") = 50 "Heat transfer, Water, Free convection";
    parameter Real hw_forced(unit = "W/m²K") = 1000 "Heat transfer, Water, Forced convection";
    parameter Real d1k1(unit = "m²K/W") = 1 / 0.5 "Typical value for 5cm insulator";
    parameter Real P0(unit = "W") = 15e3 "Maximum power, Heating element";
    // Initial Conditions
    parameter Real T0(unit = "°C") = 25 "Initial temperature";
    parameter Real H0(unit = "J/kg") = cp * (T0 - Tref) "Intial specific enthalpy";
    parameter Real U0(unit = "J") = m * H0 "Initial internal energy";
    // State variables
    Real U(start = U0, fixed = true, unit = "J") "Internal energy";
    // Auxiliary Variables
    Real T(unit = "°C") "Temperature in tank";
    Real Hm(unit = "J/kg") "Specific enthalpy in tank";
    Real Hdi(unit = "W") "Influent enthalpy flow";
    Real Hde(unit = "W") "Efluent enthalpy flow";
    Real Hmi(unit = "J/kg") "Specific enthalpy in influent flow";
    Real Qd(unit = "W") "Total heat entering the system";
    Real Qdel(unit = "W") "Heat energy from electric heater";
    Real Qda(unit = "W") "Heat energy from surroundings";
    Real md(unit = "kg/s") "Mass flow through tank";
    Real Uloss(unit = "W/m²K") "Thermal heat transfer unit to surroundings";
    Real hw(unit = "W/m²K") "Heat transfer, Water";
    Real a11;
    Real b11;
    Real b12;
    Real b13;
    Real b14;
    // Inputs
    input Real Ti;
    input Real Ta;
    input Real Vd;
    input Real up;
    input Real uv;
    // Outputs
  equation
// Differential equations
    der(U) = Hdi - Hde + Qd;
// Algebraic equations
    U = m * Hm;
    Hdi = md * Hmi;
    Hde = md * Hm;
    Hmi = cp * (Ti - Tref);
    Hm = cp * (T - Tref);
    Qd = Qdel + Qda;
    Qdel = P0 * up;
    Qda = Uloss * As * (Ta - T);
    Uloss = 1 / (1 / ha + d1k1 + 1 / hw);
    hw = if uv < 1e-7 then hw_free else hw_forced;
    md = rho * (Vd / (l2m3 * m2s)) * uv;
    a11 = -(4*Vd*uv)/(PI*d^2*h)-(4*Uloss*(d/2)+h)/(cp*rho*d*h);
    b11 = (4*Vd*uv)/(PI*d^2*h);
    b12 = (4*uv*(Ti-T))/(PI*d*h);
    b13 = (4*Vd*(Ti-T))/(PI*d^2*h);
    b14 = (4*P0)/(cp*rho*PI*d^2*h);
    
  end ModSimpleWaterHeaterMatrix;


  model ModSlicedWaterHeater
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 2
      Date: 28.10.19
      Purpose: Model of a waterheater sliced into N parts each beeing perfectly mixed
      */
    // Constants
    constant Real PI = 3.14;
    constant Real g(unit = "m/s²") = 9.81 "Gravitational contant";
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Integer N = 10 "Number of slices";
    constant Integer HeaterPosition = 7 "The heater slice position";
    // System parameters
    parameter Real d(unit = "m") = 0.5 "Tank diameter";
    parameter Real h(unit = "m") = 1.5 "Tank height";
    parameter Real A(unit = "m²") = 0.25 * PI * d ^ 2 "Tank bottom and top area";
    parameter Real dz = h / N "Slice height";
    parameter Real V(unit = "m³") = A * dz "Slice volume";
    parameter Real m(unit = "kg") = rho * V "Slice mass";
    // Thermodynamic parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    parameter Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    parameter Real ha(unit = "W/m²K") = 3 "Heat transfer, Air";
    parameter Real hw_free(unit = "W/m²K") = 50 "Heat transfer, Water, Free convection";
    parameter Real hw_forced(unit = "W/m²K") = 1000 "Heat transfer, Water, Forced convection";
    parameter Real d1k1(unit = "m²K/W") = 1 / 0.5 "Typical value for 5cm insulator";
    parameter Real kt(unit = "W/mK") = 0.6 "Thermal conductivity, Water";
    parameter Real P0(unit = "W") = 15e3 "Maximum power, Heating element";
    parameter Real kappa = 0.41 "von Karman constant";
    parameter Real alfa(unit = "1/K") = 303e-6 "Thermal expansion coefficient";
    // Initial Conditions
    parameter Real T0(unit = "°C") = 25 "Initial temperature in each slice";
    parameter Real Hm0(unit = "J/kg") = cp * (T0 - Tref) "Initial enthalpy in each slice";
    parameter Real U0(unit = "J") = m * Hm0 "Initial internal energy in each slice";
    // State variables
    Real U[N](each start = U0, each fixed = true, each unit = "J") "Internal energy in each slice";
    // Auxiliary Variables
    Real T[N](each start = T0, each unit = "°C") "Temperature in each slice";
    Real Hm[N + 1](each unit = "J/kg") "Specific entalpy in each slice including influent flow";
    Real Hd[N + 1](each unit = "W") "Entalpy flow from each slice including influent flow";
    Real Q[N](each unit = "W") "Heat entering each slice";
    Real Qel[N](each unit = "W") "Heat from heater";
    Real Qa[N](each unit = "W") "Heat from surroundings";
    Real Qd[N](each unit = "W") "Heat diffusion for each slice";
    Real Qdd[N + 1](each unit = "W/m²") "Heat diffusion flux for each slice";
    Real Qdb[N + 1](each unit = "W/m²") "Boyant turbulence mixing for each slice";
    Real dTdz[N + 1](each unit = "°C/m") "Temperature gradient in tank";
    Real d2Tdz2[N + 1](each unit = "°C/m²") "Temperature hessian in tank";
    Real dQdbdz[N + 1](each unit = "W/m³") "Boyant turbulence mixing gradient";
    Real As[N](each unit = "m²") "Tank wall surface area for each slice";
    Real md(unit = "kg/s") "Total mass flow through tank";
    Real Uloss(unit = "W/m²K") "Thermal heat transfer unit to surroundings";
    Real hw(unit = "W/m²K") "Heat transfer, Water";
    Real kb[N + 1] "Boyant conductivity, Water";
    Real Pk[N] "Position array for heater placement";
    // Debug variables
    Real sumU, sumderU, sumHi, sumHe, sumQ, sumQa, sumQel, deltaU;
    // Inputs
    input Real Ti(unit = "°C") "Influent temperatur";
    input Real Ta(unit = "°C") "Ambient temperature";
    input Real Vd(unit = "l/min") "Total volum flow";
    input Real up "Output to heater, 0.0-1.0";
    input Real uv "Output to shunt valve, 0.0-1.0";
    // Outputs
    output Real Te;
  equation
// Differential equations
    der(U[:]) = Hd[1:N] - Hd[2:N + 1] + Q[:];
// Algebraic equations
    U[:] = m * Hm[2:N + 1];
    Hd[:] = md * Hm[:];
    Hm[2:N + 1] = cp * (T[:] .- Tref);
    Hm[1] = cp * (Ti - Tref);
    Q[:] = Qa[:] + Qel[:] + Qd[:];
    Qel[:] = P0 * up * Pk[:];
    Qa[:] = Uloss .* As[:] .* (Ta .- T[:]);
    Qd[:] = (Qdd[1:N] + Qdb[1:N] - Qdd[2:N + 1] - Qdb[2:N + 1]) * A;
    Qdd[:] = -kt * dTdz[:];
    Qdb[N + 1] = 0;
    Qdb[1] = 0;
    dQdbdz[2:N] = -kb[2:N] .* d2Tdz2[2:N];
    As[1] = PI * d * dz + A;
    As[2:N - 1] = PI * d * dz * ones(N - 2);
    As[N] = PI * d * dz + A;
    Uloss = 1 / (1 / ha + d1k1 + 1 / hw);
    hw = if uv < 1e-7 then hw_free else hw_forced;
    md = rho * (Vd / (l2m3 * m2s)) * uv;
    for i in 1:N + 1 loop
      kb[i] = if dTdz[i] < 0.0 then rho * cp * (kappa * d) ^ 2 * sqrt(g * alfa * abs(dTdz[i])) else 0.0;
    end for;
    Te = T[N];
// Finite difference approximation
    dTdz[1] = 0.0;
    dTdz[2:N] = (T[2:N] - T[1:N - 1]) / dz;
    dTdz[N + 1] = 0.0;
    d2Tdz2[1] = 0.0;
    d2Tdz2[2:N] = (dTdz[2:N] - dTdz[1:N - 1]) / dz;
    d2Tdz2[N + 1] = 0.0;
    dQdbdz[1] = 0.0;
    dQdbdz[2:N] = (Qdb[2:N] - Qdb[1:N - 1]) / dz;
    dQdbdz[N + 1] = 0.0;
// Position functions
    Pk[1:HeaterPosition - 1] = zeros(HeaterPosition - 1);
    Pk[HeaterPosition + 1:N] = zeros(N - HeaterPosition);
    Pk[HeaterPosition] = 1.0;
// Sum functions (debugging)
    sumU = U[:] * ones(N);
    sumderU = der(U[:]) * ones(N);
    sumHi = Hd[1];
    sumHe = Hd[N + 1];
    sumQ = Q[:] * ones(N);
    sumQa = Qa[:] * ones(N);
    sumQel = Qel[:] * ones(N);
    deltaU = sumQ + sumHi - sumHe;
  end ModSlicedWaterHeater;

  model ModSlicedWaterHeater2
    /*
      Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
      Revision: 1
      Date: 28.10.19
      Purpose: Model of a waterheater sliced into N parts each beeing perfectly mixed
      */
    // Constants
    constant Real PI = 3.14;
    constant Real g(unit = "m/s²") = 9.81 "Gravitational contant";
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Integer N = 20 "Number of slices";
    constant Real hPosition = 1.15 "The heater position";
    // System parameters
    parameter Real d(unit = "m") = 0.5 "Tank diameter";
    parameter Real h(unit = "m") = 1.5 "Tank height";
    parameter Real A(unit = "m²") = 0.25 * PI * d ^ 2 "Tank bottom and top area";
    parameter Real dz = h / N "Slice height";
    parameter Real V(unit = "m³") = A * dz "Slice volume";
    parameter Real m(unit = "kg") = rho * V "Slice mass";
    parameter Real hHeight = dz * 0.01 "The heater position";
    // Thermodynamic parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    parameter Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    parameter Real ha(unit = "W/m²K") = 3 "Heat transfer, Air";
    parameter Real hw_free(unit = "W/m²K") = 50 "Heat transfer, Water, Free convection";
    parameter Real hw_forced(unit = "W/m²K") = 1000 "Heat transfer, Water, Forced convection";
    parameter Real d1k1(unit = "m²K/W") = 1 / 0.5 "Typical value for 5cm insulator";
    parameter Real kt(unit = "W/mK") = 0.6 "Thermal conductivity, Water";
    parameter Real P0(unit = "W") = 15e3 "Maximum power, Heating element";
    parameter Real kappa = 0.41 "von Karman constant";
    parameter Real alfa(unit = "1/K") = 303e-6 "Thermal expansion coefficient";
    // Initial Conditions
    parameter Real T0(unit = "°C") = 25 "Initial temperature in each slice";
    parameter Real Hm0(unit = "J/kg") = cp * (T0 - Tref) "Initial enthalpy in each slice";
    parameter Real U0(unit = "J") = m * Hm0 "Initial internal energy in each slice";
    // State variables
    Real U[N](each start = U0, each fixed = true, each unit = "J") "Internal energy in each slice";
    // Auxiliary Variables
    Real T[N](each start = T0, each unit = "°C") "Temperature in each slice";
    Real Hm[N + 1](each unit = "J/kg") "Specific entalpy in each slice including influent flow";
    Real Hd[N + 1](each unit = "W") "Entalpy flow from each slice including influent flow";
    Real Q[N](each unit = "W") "Heat entering each slice";
    Real Qel[N](each unit = "W") "Heat from heater";
    Real Qa[N](each unit = "W") "Heat from surroundings";
    Real Qd[N](each unit = "W") "Heat diffusion for each slice";
    Real Qdd[N - 1](each unit = "W/m²") "Heat diffusion flux for each slice";
    //Real _Qdb[N-1](each unit="W/m²")         "Boyant turbulence mixing for each slice";
    Real Qdb[N - 1](each unit = "W/m²") "Boyant turbulence mixing for each slice";
    Real dTdz[N - 1](each unit = "°C/m") "Temperature gradient in tank";
    Real _dTdz[N](each unit = "°C/m") "Temperature gradient in tank";
    Real d2Tdz2[N](each unit = "°C/m²") "Temperature hessian in tank";
    //Real _d2Tdz2[N-2](each unit="°C/m²")     "Temperature hessian in tank";
    Real dQdbdz[N](each unit = "(W/m²)/m") "Boyant turbulence mixing gradient";
    //Real _dQdbdz[N-1](each unit="(W/m²)/m")      "Boyant turbulence mixing gradient";
    Real As[N](each unit = "m²") "Tank wall surface area for each slice";
    Real md(unit = "kg/s") "Total mass flow through tank";
    Real Uloss(unit = "W/m²K") "Thermal heat transfer unit to surroundings";
    Real hw(unit = "W/m²K") "Heat transfer, Water";
    Real kb[N] "Boyant conductivity, Water";
    Real Pk[N] "Position array for heater placement";
    // Debug variables
    Real sumU, sumderU, sumHi, sumHe, sumQ, sumQa, sumQel, deltaU;
    // Inputs
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 4 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 13 "Total volum flow";
    parameter Real up = 1 "Output to heater, 0.0-1.0";
    parameter Real uv = 0 "Output to shunt valve, 0.0-1.0";
    // Outputs
    output Real Te;
  equation
// Differential equations
    der(U[:]) = Hd[1:N] - Hd[2:N + 1] + Q[:];
// Algebraic equations
    U[:] = m * Hm[2:N + 1];
    Hd[:] = md * Hm[:];
    Hm[2:N + 1] = cp * (T[:] .- Tref);
    Hm[1] = cp * (Ti - Tref);
    Q[:] = Qa[:] + Qel[:] + Qd[:];
    Qel[:] = up * Pk[:];
    Qa[:] = Uloss .* As[:] .* (Ta .- T[:]);
//Qd[1] = -Qdd[1]*A;
//Qd[2:N-1] = (Qdd[1:N-2] - Qdd[2:N-1])*A;
//Qd[N] =  Qdd[N-1]*A;
    Qd[1] = (Qdd[1] + Qdb[1]) * A;
    Qd[2:N - 1] = (Qdd[1:N - 2] + Qdb[1:N - 2] - Qdd[2:N - 1] - Qdb[2:N - 1]) * A;
    Qd[N] = ((-Qdd[N - 1]) - Qdb[N - 1]) * A;
    Qdd[:] = -kt * dTdz[:];
//Qdd[N] = 0.0;
//dQdbdz[N] = 0;
    dQdbdz[:] = -kb[:] .* d2Tdz2[:];
//dQdbdz[N] = 0;
//dQdbdz[1:N-1] = -kb[1:N-1].*d2Tdz2[1:N-1];
    As[1] = PI * d * dz + A;
    As[2:N - 1] = PI * d * dz * ones(N - 2);
    As[N] = PI * d * dz + A;
    Uloss = 1 / (1 / ha + d1k1 + 1 / hw);
    hw = if uv < 1e-7 then hw_free else hw_forced;
    md = rho * (Vd / (l2m3 * m2s)) * uv;
    for i in 1:N loop
      kb[i] = if _dTdz[i] < 0.0 then rho * cp * (kappa * d) ^ 2 * sqrt(g * alfa * abs(_dTdz[i])) else 0.0;
    end for;
    Te = T[N];
// Finite difference approximation
    dTdz[:] = (T[2:N] - T[1:N - 1]) / dz;
    _dTdz[1] = dTdz[1];
    _dTdz[2:N - 1] = (dTdz[1:N - 2] + dTdz[2:N - 1]) / 2;
    _dTdz[N] = dTdz[N - 1];
//_d2Tdz2[:] = (dTdz[2:N-1]-dTdz[1:N-2])/(dz);
    d2Tdz2[1] = dTdz[1] / dz;
    d2Tdz2[2:N - 1] = (dTdz[2:N - 1] - dTdz[1:N - 2]) / dz;
    d2Tdz2[N] = -dTdz[N - 1] / dz;
    Qdb[1] = dQdbdz[1] * dz;
    Qdb[2:N - 1] = Qdb[1:N - 2] + dQdbdz[2:N - 1] * dz;
//_Qdb[1] = 0;
//dQdbdz[1] = _Qdb[1]/(dz);
//dQdbdz[2:N-1] = (_Qdb[2:N-1]-_Qdb[1:N-2])/(dz);
//dQdbdz[N] = -Qdb[N-1]/(dz);
//_dQdbdz[N-1] = (Qdb[2:N-1]-Qdb[1:N-2])/(dz);
//dQdbdz[N] = (Qdb[N+1]-Qdb[N])/(dz);
// Position functions
    for i in 1:N loop
      if (i - 1) * dz < hPosition and i * dz > hPosition then
        Pk[i] = min(P0 * ((i * dz - hPosition) / hHeight), P0);
      elseif (i - 1) * dz > hPosition and i * dz < hPosition + hHeight then
        Pk[i] = P0 * (dz / hHeight);
      elseif (i - 1) * dz < hPosition + hHeight and i * dz > hPosition + hHeight then
        Pk[i] = P0 * ((hPosition + hHeight - (i - 1) * dz) / hHeight);
      else
        Pk[i] = 0;
      end if;
    end for;
// Sum functions (debugging)
    sumU = U[:] * ones(N);
    sumderU = der(U[:]) * ones(N);
    sumHi = Hd[1];
    sumHe = Hd[N + 1];
    sumQ = Q[:] * ones(N);
    sumQa = Qa[:] * ones(N);
    sumQel = Qel[:] * ones(N);
    deltaU = sumQ + sumHi - sumHe;
  end ModSlicedWaterHeater2;

  model ModSimpleHeater
    /*
          Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
          Revision: 1
          Date: 25.10.19
          Purpose: Model of a perfectly mixed water heater
          */
    // Constants
    constant Real PI = 3.14;
    constant Real l2m3 = 1000 "Conversion from liter to m³";
    constant Real m2s = 60 "Conversion from minutes to seconds";
    constant Real cp(unit = "kJ/kgK") = 4190 "Heat capacity of water";
    // Parameters
    parameter Real Tref(unit = "°C") = 25 "Reference temperature";
    parameter Real rho(unit = "kg/m³") = 1000 "Water density";
    parameter Real d(unit = "m") = 0.25 "Tank diameter";
    parameter Real h(unit = "m") = 0.5 "Tank height";
    parameter Real A(unit = "m²") = 0.25 * PI * d ^ 2 "Tank bottom and top area";
    parameter Real V(unit = "m³") = A * h "Tank volume";
    parameter Real As(unit = "m²") = 2 * A + PI * d * h "Tank surface area";
    parameter Real m(unit = "kg") = rho * V "Water mass in tank";
    parameter Real ha(unit = "W/m²K") = 9.3 "Heat transfer, Air";
    parameter Real hw_free(unit = "W/m²K") = 50 "Heat transfer, Water, Free convection";
    parameter Real hw_forced(unit = "W/m²K") = 1000 "Heat transfer, Water, Forced convection";
    parameter Real d1k1(unit = "m²K/W") = 1 / 100 "Typical value for 5cm insulator";
    parameter Real P0(unit = "W") = 2.5e3 "Maximum power, Heating element";
    // Initial Conditions
    parameter Real T0(unit = "°C") = 15 "Initial temperature";
    parameter Real H0(unit = "J/kg") = cp * (T0 - Tref) "Intial specific enthalpy";
    parameter Real U0(unit = "J") = m * H0 "Initial internal energy";
    
    parameter Real Hevap(unit = "J/kg") = 800  "Latent heat of evaporization";
    parameter Real vair(unit = "m/s") = 0 "air speed";
    parameter Real xamb(unit = "kg/kg") = 10.41e-3 "specific humidity";
    parameter Real Pa(unit = "Pa") = 1.0135e5 "ambient pressure";
    
    // State variables
    Real U(start = U0, fixed = true, unit = "J") "Internal energy";
    // Auxiliary Variables
    
    Real Ti(unit = "°C") "Influent temperatur";
    Real T(unit = "°C") "Temperature in tank";
    Real Hm(unit = "J/kg") "Specific enthalpy in tank";
    Real Hdi(unit = "W") "Influent enthalpy flow";
    Real Hde(unit = "W") "Efluent enthalpy flow";
    Real Hmi(unit = "J/kg") "Specific enthalpy in influent flow";
    Real Qd(unit = "W") "Total heat entering the system";
    Real Qdel(unit = "W") "Heat energy from electric heater";
    Real Qda(unit = "W") "Heat energy from surroundings";
    Real Qdevap(unit = "W") "Heat energy from evaporization";
    Real xsat(unit = "kg/kg") "Specific saturation humidity";
    Real Psat(unit = "Pa") "vapour saturation pressure";
    Real md(unit = "kg/s") "Mass flow through tank";
    Real Uloss(unit = "W/m²K") "Thermal heat transfer unit to surroundings";
    Real hw(unit = "W/m²K") "Heat transfer, Water";
    // Inputs
    parameter Real Ta(unit = "°C") = 15 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 10 "Total volum flow";
    input Real up "Output to heater, 0.0-1.0";
    parameter Real uv = 1 "Output to shunt valve, 0.0-1.0";
    // Outputs
    //output Real T;
  equation
// Differential equations
    der(U) = Hdi - Hde + Qd;
// Algebraic equations
    U = m * Hm;
    Hdi = md * Hmi;
    Hde = md * Hm;
    Hmi = cp * (Ti - Tref);
    Hm = cp * (T - Tref);
    Qd = Qdel + Qda - Qdevap;
    Qdel = P0 * up;
    Qda = Uloss * As * (Ta - T);
    Uloss = 1 / (1 / ha + d1k1 + 1 / hw);
    hw = if uv < 1e-7 then hw_free else hw_forced;
    md = rho * (Vd / (l2m3 * m2s)) * uv;
    Ti = T;
    
    Qdevap = (25 + 19*vair)*A*(xsat-xamb)*Hevap;
    xsat = (0.663*Psat)/(Pa-0.378*Psat);
    Psat = (6.1078*exp((17.2693882*((T+273.16)-273.16))/((T+273.16)-35.86)))*100;
    
  end ModSimpleHeater;

  model SimSimpleHeater
    /*
          Authors: Kristian Dyb Strand, Kristian Aasbø Hansen and Erik Kristoffer Rummelhoff
          Revision: 1
          Date: 25.10.19
          Purpose: Simulation of model "ModSimpleWaterHeater"
          */
    // Instatiating model
    ModSimpleHeater swh;
    Modelica.Blocks.Continuous.LimPID PI(
      k = 0.732,
      Ti=240,
      yMax = 1,
      Ni=0.1,
      limitsAtInit=false,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Td=0.1);
    // Parameters
    parameter Real Ti(unit = "°C") = 27 "Influent temperatur";
    parameter Real Ta(unit = "°C") = 15 "Ambient temperature";
    parameter Real Vd(unit = "l/min") = 10 "Total volum flow";
    Real up "Output to heater, 0.0-1.0";
    parameter Real uv = 1;
    Real sp;
    parameter Real kp = 15;
    // Input variables
    //Real uv "Output to shunt valve, 0.0-1.0";
    // Output variable
    Real Te(unit = "°C") "Effluent temperature";
  equation
    //swh.Ti = Te;
    //swh.Ta = Ta;
    //swh.Vd = Vd;
//uv = if time < 3600 then 0 else 1;
    //swh.uv = uv;
    PI.u_m = Te;
    PI.u_s = sp;
    up = PI.y;
    sp = if time <3600 then 60 else 63;
    //up = if time <3600 then 0.23 else 0.25;
    swh.up = up;
    Te = swh.T;
  end SimSimpleHeater;
end projectWaterHeater;
