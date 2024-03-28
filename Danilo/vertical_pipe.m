clear
close
clc

tic
global  G dL D tk r1 r2 A eps ks ts;
% Pipe parameters
D = 0.1;            % [m] pipe diameter
tk = 0.005;         % [m] pipe shell thickness
r1 = D/2;           % [m] pipe internal radius
r2 = r1 + tk;       % [m] pipe outer radius
A = pi/4*D^2;       % [m^2] pipe cross section area 
eps = 2*10^(-5);    % [m] pipe wall roughness
ks = 20;            % [W/(mK)] pipe solid shell conductivity

% Boundary conditions
ts = 200.0;         % [°C] pipe shell outer temperature

% Inflow conditions
Re_in = 10^5;         % inflow Reynolds number
t_in = 20.0;          % [°C] inflow temperature
p_in = 10^5;          % [Pa] inflow pressure

% derived flow quantities                   
i_sat = XSteam('hL_p',p_in/10^5)*1000;                                          % [J/kg] water saturation enthaply at inlet pressure
i_in = XSteam('h_pT',p_in/10^5,t_in)*1000;                                           % [J/kg] water inlet enthalpy
mu_in = XSteam('my_pT',p_in/10^5,t_in);                                      % [Pa s] water inlet dynamic viscosity 
rho_in = XSteam('rho_pT',p_in/10^5,t_in);                                             % [kg/m^3] water inlet density

G = A*mu_in/D*Re_in;                                                    % [kg/s] mass flow rate

% Finite volume length
dL = 0.01;       % [m]


% Solution initialization
i = [i_in i_in];            % [J/kg] water enthalpy
tf = [t_in t_in];           % [°C] water temperature
x = 0;                 % [m] pipe overall length
p = [p_in p_in];            % [Pa] pressure

% Iteration loop:
% at each step, it calculates the temperature and the enthalpy at that temperature: 
% if it is still lower than saturation enthalpy at that pressure it continues to iterate
% adding a finite volume to the pipe
while i(end) < i_sat && p(end) > 10^3

    U = thermal_resistance(p(end),tf(end));
    % find enthalpy at the next step from balance equation
    i1 = i(end) + U/G*(ts-tf(end));
    % find pressure drop
    p1 = p(end) - pressure_drop(p(end),tf(end)) - XSteam('rho_pT',p(end)/10^5,tf(end))*9.81*dL;
    % find temperature at the next step
    t1 = XSteam('T_ph',p(end)/10^5,i(end)/1000);
    % update saturation enthalpy due to pressure drop
    i_sat = XSteam('hL_p',p(end)/10^5)*1000;

    % save values
    i = [i, i1];
    tf = [tf, t1];
    x = x+dL;
    p = [p, p1];
end

toc

% function that calculates pressure drop using water conditions at given pressure and temperature
function dp = pressure_drop(p,t)
    global  G dL D tk r1 r2 A eps ks ts;
    % extract water properties
    rho = XSteam('rho_pT',p/10^5,t);            % [kg/m^3] density
    mu = XSteam('my_pT',p/10^5,t);              % [Pa s] dynamic viscosity

    v = G/(rho*A);                   % [m/s] average bulk velocity (mass flow is retained)
    Re = D*rho*v/mu;                 % Reynolds number for these water conditions

    % calculate pressure drop using Colebrook correlation
    f = (10^-3:10^-3:10);                                                % range of f values
    left_term = 1./(f.^0.5);                                             % left term of Colebrook correlation
    right_term = -2*log(eps/(3.7*D)+2.51./(Re.*f.^0.5));                 % right term of Colebrook correlation
    [~, min_index] = min(abs(left_term-right_term));                     % find the index where the difference between the two terms is near 0

    fD = 1./(left_term(min_index).^2)*4;                                 % Darcy-Weisbach coefficient calculation
    dp = dL*fD*rho/2*v^2/D;                                             % Darcy-Weisbach pressure drop
end

% overall thermal resistance
function u = thermal_resistance(p,t)
    global  G dL D tk r1 r2 A eps ks ts;
    % water properties at found conditions
    rho = XSteam('rhoL_T',t);            % [kg/m^3] density
    mu = XSteam('my_pT',p/10^5,t);              % [Pa s] dynamic viscosity
    cp = XSteam('Cp_pT',p/10^5,t)*1000;             % [J/(kg K)] specific heat capacity
    k = XSteam('tc_pT',p/10^5,t);                    % [W/(m K)] water conductivity

    v = G/(rho*A);                       % [m/s] average bulk velocity (mass flow is retained)
    Re = D*rho*v/mu;                     % Reynolds number for these water conditions
    Pr = mu*cp/k;                        % Prandtl number for these water conditions
    Nu = 0.023*Re^0.8*Pr^0.4;            % Nusslet number for these water conditions (Dittus–Boelter)
    h = Nu*k/D;                          % [W/(m^2 K)] convective heat transfer coefficient

    u = (1./(h*2*pi*r1*dL)+log(r2/r1)./(2*pi*ks*dL))^(-1);    % [W/K] overall thermal resistance
end
