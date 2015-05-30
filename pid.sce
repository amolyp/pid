// Circuit cellar #221 - The darker side : PID controls
// Illustration scripts to be executed under Scilab V4.1.2
// Robert Lacoste

//---------------------------------------------------------------------------
// Physical parameters
//---------------------------------------------------------------------------

// Hypothesis : 
// - The heater has an homogeneous temperature, and heats the tube
// - The tube has an homogeneous temperature, heats the tip and dissipates
// - The tip has an homogeneous temperature, and dissipates.

// Ambiant
T_ambiant=20;               // ambiant temperature (°C)

// Heater
m_heater=20;                // mass of the heater (g)
c_heater=0.385;             // Specific heat capacity of copper (J per g per °C)
C_heater=m_heater*c_heater; // Heat capacity of the heater (J per °C)
p_heat=80;                  // Heater power (W)
k_heater_to_tube=1.0;       // Heat transfer coefficient from heater to tube (W per °C)

// Tube
m_tube=20;                  // mass of the tube (g)
c_tube=0.385;               // Specific heat capacity of copper (J per g per °C)
C_tube=m_tube*c_tube;       // Heat capacity of the tube (J per °C)
k_tube_to_tip=0.4;          // Heat transfer coefficient from tube to tip (W per °C)
kloss_tube=0.04;            // Heat loss coefficient from tube to ambiant (W per °C)

// Tip
m_tip=10;                   // mass of the tip (g)
c_tip=0.385;                // Specific heat capacity of copper (J per g per °C)
C_tip=m_tip*c_tip;          // Heat capacity of the tip (J per °C)
kloss_tip=0.02;             // Heat loss coefficient from tip to ambiant (W per °C)

//---------------------------------------------------------------------------
// Simulation parameters
//---------------------------------------------------------------------------

nsteps=3000;                // number of simulation steps
dt=0.2;                     // time step for simulation (s)
t=0:dt:nsteps*dt;           // Time vector
T_heater=0*t;               // Create vector of heater temperatures
T_tube=0*t;                 // Create vector of tube temperatures
T_tip=0*t;                  // Create vector of tip temperatures
Command=0*t;                // Create vector of commands (0 to 1)
T_target=350;               // Target tip temperature (°C)

//---------------------------------------------------------------------------
// Simulation
//---------------------------------------------------------------------------

T_heater(1)=T_ambiant;                                // Initial conditions
T_tube(1)=T_ambiant;
T_tip(1)=T_ambiant;
err=0;
integral_err=0;

for i=2:length(t)                                     // Main loop through time

  // ===================
  // Command calculation
  //====================
  
  kp=0.1;
  kd=1;
  ki=0.0005;
  maxintegral=1000;
  previous_err=err;
  integral_err=integral_err+err;
  integral_err=max(min(integral_err,maxintegral),-maxintegral);
  err=T_target-T_tip(i-1);
  c = kp*err + kd*(err-previous_err) + ki*integral_err;
  c=max(min(c,1),0);
  Command(i)=c;
    
  // ===================
  // Thermal simulation
  //====================
  
  Pgain=p_heat*Command(i-1);                                      // Heater gains power if powered
  Pheater_to_tube = k_heater_to_tube*(T_heater(i-1)-T_tube(i-1)); // Heat transfer to tube
  Ptube_to_tip = k_tube_to_tip*(T_tube(i-1)-T_tip(i-1));          // Heat transfer to tip
  Ploss_tube = kloss_tube*(T_tube(i-1)-T_ambiant);                // Heat loss from tube to ambiant
  Ploss_tip = kloss_tip*(T_tip(i-1)-T_ambiant);                   // Heat loss from tip to ambiant
  deltaT_heater=(Pgain-Pheater_to_tube)*dt/C_heater;
  T_heater(i)=T_heater(i-1)+deltaT_heater;                        // Heater temperature change
  deltaT_tube=(Pheater_to_tube-Ptube_to_tip-Ploss_tube)*dt/C_tube;
  T_tube(i)=T_tube(i-1)+deltaT_tube;                              // Tube temperature change
  deltaT_tip=(Ptube_to_tip-Ploss_tip)*dt/C_tip;
  T_tip(i)=T_tip(i-1)+deltaT_tip;                                 // Tip temperature change
  
end;

// ===================
// Plot results
// ===================

subplot(3,1,1); 
xtitle('PID control; Kp=0.1, Kd=1, Ki=0.0005');
plot2d(t,T_heater,rect=[1,0,t(nsteps),max(T_heater)],style=[6]); 
plot2d(t,T_tube,rect=[1,0,t(nsteps),max(T_heater)],style=[3]); 
plot2d(t,T_tip,rect=[1,0,t(nsteps),max(T_heater)],style=[2]); 
subplot(3,1,2); 
plot2d(t,Command,rect=[1,-0.1,t(nsteps),1.1],style=[5]);
subplot(4,1,4); 
plot2d(t(2000:3000),T_tip(2000:3000),rect=[t(2000),min(T_tip(2000:3000)),t(3000),max(T_tip(2000:3000))],style=[2]);
halt();
xdel();

