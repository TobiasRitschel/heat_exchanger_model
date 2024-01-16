% Simulate heat exchanger (relevant to a nuclear reactor with circulated fuel)
% Clear command window
clc;

% Clear variables, etc.
clear all;

% Close figures, etc.
close all;

% Remove added paths
restoredefaultpath;

%% System
% Right-hand side function
f = @heat_exchanger_right_hand_side;

% Spatial domain
z0 = 0;
zf = 1;

% Temporal domain
t0 =  0;
tf = 10;

% Collect system parameters
p.z0 = z0;
p.zf = zf;

%% Parameters
% Radii
p.Rf = 0.49;
p.Rt = 0.51;
p.Rc = 1.00;

% Heat capacities
p.cPf = 1;
p.cPt = 1;
p.cPc = 1;

% Heat conductivities
p.ktf = 1;
p.ktc = 1;

% Concentrations (identical to densities for this simple model)
p.cf = 1;
p.ct = 1;
p.cc = 1;

% Standard temperature (no effect on simulation)
p.T0 = 273.15; % [K]

% Standard enthalpies (no effect on simulation)
p.hf0 = 0;
p.hc0 = 0;

%% Manipulated inputs
% Flow velocities
vf =  1;
vc = -1;

% Collect manipulated inputs
u = [vf; vc];

%% Disturbance variables
% Inlet temperatures
Tfin = 700; % [K]
Tcin = 600; % [K]

% Collect disturbance variables
d = [Tfin; Tcin];

%% Discretization parameters and options
% Number of time steps (only for plotting)
N = 1000;

% Number of finite volumes
M = 100;

% Collect discretization parameters
p.M = M;

% Time span
tspan = linspace(t0, tf, N+1);

% Options
fsolve_opts = optimoptions('fsolve',    ...
    'FunctionTolerance',        1e-12,  ...
    'Display',                  'None');

% Options
ode_opts = [];

%% Steady state
% Initial guess of steady state
xinit(    1:  M, 1) =      Tfin;
xinit(  M+1:2*M, 1) = 0.5*(Tfin + Tcin);
xinit(2*M+1:3*M, 1) =             Tcin;

% Time (not relevant)
t = [];

% Solve for the steady state
x0 = fsolve(@steady_state_wrapper, xinit, fsolve_opts, ...
    f, t, u, d, p);

%% Visualize steady state solution
% Create figure
figure(1);

% Visualize steady state
visualize_steady_state(x0, p);

%% Simulation (increased fuel inlet temperature)
% Increase fuel inlet temperature
Tfin = 800; % [K]

% Update vector with disturbance variables
d(1) = Tfin;

% Simulate
[t, X] = ode15s(f, tspan, x0, ode_opts, ...
    u, d, p);

%% Visualize simulation
% Create figure
figure(2);

% Visualize simulation
visualize_simulation(t, X, p);

%% Simulation (Decreased coolant inlet temperature)
% Reset fuel inlet temperature
Tfin = 700; % [K]

% Decrease coolant inlet temperature
Tcin = 500; % [K]

% Update vector with disturbance variables
d(1) = Tfin;
d(2) = Tcin;

% Simulate
[t, X] = ode15s(f, tspan, x0, ode_opts, ...
    u, d, p);

%% Visualize simulation
% Create figure
figure(3);

% Visualize simulation
visualize_simulation(t, X, p);

%% Simulation (Increased fuel velocity)
% Update fuel velocity
vf = 5;

% Reset coolant inlet temperature
Tcin = 600; % [K]

% Update vector with manipulated inputs
u(1) = vf;

% Update vector with disturbance variables
d(2) = Tcin;

% Simulate
[t, X] = ode15s(f, tspan, x0, ode_opts, ...
    u, d, p);

%% Visualize simulation
% Create figure
figure(4);

% Visualize simulation
visualize_simulation(t, X, p);

%% Simulation (Increased fuel velocity)
% Reset fuel velocity and update coolant velocity
vf =  1;
vc = -5;

% Update vector with manipulated inputs
u(1) = vf;
u(2) = vc;

% Simulate
[t, X] = ode15s(f, tspan, x0, ode_opts, ...
    u, d, p);

%% Visualize simulation
% Create figure
figure(5);

% Visualize simulation
visualize_simulation(t, X, p);

%% Functions
function visualize_steady_state(x0, p)
% Extract parameters
z0 = p.z0;
zf = p.zf;

M = p.M;

% Spatial coordinates
z = linspace(z0, zf, M+1);

% Extract temperatures
Tfs = x0(      [1:M, M]);
Tts = x0(  M + [1:M, M]);
Tcs = x0(2*M + [1:M, M]);

% Visualize steady state fuel temperatures
stairs(z, Tfs);

% Keep adding to the plot
hold on;

% Visualize steady state tubing temperatures
stairs(z, Tts);

% Visualize steady state coolant temperatures
stairs(z, Tcs);

% Stop adding to the plot
hold off;

% Axis labels
xlabel('Spatial coordinate, z');
ylabel('Temperature, T');
end

function visualize_simulation(t, X, p)
% Extract parameters
z0 = p.z0;
zf = p.zf;

M = p.M;

% Spatial coordinates
z = linspace(z0, zf, M+1);

% Extract temperatures
Tf = X(:,       [1:M, M]);
Tt = X(:,   M + [1:M, M]);
Tc = X(:, 2*M + [1:M, M]);

% Visualize steady state fuel temperatures
surf(z, t, Tf);

% Keep adding to the plot
hold on;

% Visualize steady state tubing temperatures
surf(z, t, Tt);

% Visualize steady state coolant temperatures
surf(z, t, Tc);

% Stop adding to the plot
hold off;

% Turn off grid lines
shading interp;

% Axis labels
xlabel('Spatial coordinate, z');
ylabel('Time, t');
zlabel('Temperature, T');
end

function rhs = steady_state_wrapper(x, f, t, u, d, p)
rhs = f(t, x, u, d, p);
end

function f = heat_exchanger_right_hand_side(~, x, u, d, p)
% Extract parameters
Rf  = p.Rf;
Rt  = p.Rt;
Rc  = p.Rc;

z0  = p.z0;
zf  = p.zf;

cPf = p.cPf;
cPt = p.cPt;
cPc = p.cPc;

ktf = p.ktf;
ktc = p.ktc;

cf  = p.cf;
ct  = p.ct;
cc  = p.cc;

T0  = p.T0;
hf0 = p.hf0;
hc0 = p.hc0;

% Number of volumes
M = p.M;

% Extract states
Tf = x(    1:  M);
Tt = x(  M+1:2*M);
Tc = x(2*M+1:3*M);

% Extract manipulated inputs
vf = u(1);
vc = u(2);

% Extract disturbance variables
Tfin = d(1);
Tcin = d(2);

% Length of finite volumes
dz = (zf - z0)/M;

% Compute areas
Af = pi*Rf^2;
At = pi*(Rt^2 - Rf^2);
Ac = pi*(Rc^2 - Rt^2);

% Compute volumes
Vf = Af*dz;
Vt = At*dz;
Vc = Ac*dz;

% Molar numbers
nf = cf*Vf;
nt = ct*Vt;
nc = cc*Vc;

% Evaluate fuel and coolant fluxes
Nf = vf*cf;
Nc = vc*cc;

% Evaluate molar enthalpies
hf = hf0 + cPf*([Tfin; Tf] - T0);
hc = hc0 + cPc*([Tc; Tcin] - T0);

% Evaluate enthalpy fluxes
Hf = Nf*hf;
Hc = Nc*hc;

% Evaluate heat transfers
Qtf = -ktf*(Tf - Tt);
Qtc = -ktc*(Tc - Tt);

Qft = -Qtf;
Qct = -Qtc;

% Evaluate right-hand sides
f(    1:  M, 1) = (-Af*(Hf(2:end) - Hf(1:end-1)) + Vf* Qtf       )/(cPf*nf);
f(  M+1:2*M, 1) = (                                Vt*(Qft + Qct))/(cPt*nt);
f(2*M+1:3*M, 1) = (-Ac*(Hc(2:end) - Hc(1:end-1)) + Vc*       Qtc )/(cPc*nc);
end