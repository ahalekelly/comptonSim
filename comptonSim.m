function y = Project %Allows M-file to be called as a function in MATLAB

%Imput Specific Design Values
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle

m_p = ****; %mass of pellet ( kg )
m_s = ****; %mass of spring piston + 1/3 spring mass ( kg )
r_piston = ****; %radius of spring piston in air gun ( m )
c = ****; %caliber ( m )
R = ****; %twist rate of rifling ( 1/m )
f_e = ****; %experimentally quantified elastic compression friction force for pellet-style ( N )
mu = ****; %kinetic coeficient of friction between pellet and barrel materials
I_p = ****; %mass moment of inertia of pellet along its lognitudinal central axis ( kg*m^2 )
g = 9.801; %force of local gravitational feild (use standard earth value for analysis) ( m/(kg*s^2) )
P_o = 1.013*10^5; %local atmospheric pressure (normally set at see level for analysis) ( Pa )
V_o = ****; %Volume of air in cylinder when the air gun is "cocked" ( m^3 )
k = ****; %spring constant ( N/m )
x_s_f = ****; %compression of spring assembly after cocking ( m )
x_s_o = ****; %compression of spring assembly before cocking ( m )
barrel_length = ****; % ( m )
R_angle = atan(R*(pi*c)); %Twist angle of Rifling Formula -- Don't Change Formula !!!

%First solve the initial value problem using "ode45"
tspan=[0,0.01];x0=[0;0]; %Time Span ( tspan ) may have to be increased or decreased for different models
[t,x] = ode45(@problem_func,tspan,x0);
x(:,1) = real((x(:,1) + ((x_s_f-x_s_o)-x(:,1)).*(x(:,1) >= (x_s_f-x_s_o))));x(:,2) = real(x(:,2));

%Make a plot for (t,x_p)
plot(t,x(:,2))
xlabel('t')
ylabel('x_p')
title('Position of Pellet in Barrel as a Function of Time')
shg

pause

%Next make a plot for (x_p,v_p)
vel_s = real((x(:,1)<(x_s_f-x_s_o)).*(k.*(x_s_f - x(:,1)) - ( (pi.*r_piston^2.*P_o.*V_o)
./(V_o - pi.*r_piston.^2.*x(:,1) + (pi./4).*c.^2.*x(:,2)) ) + pi.*r_piston.^2.*P_o).*((t)./(m_s)));
vel_p = real(sqrt( (SpringEnergy(x(:,1)) - BoreFriction(x(:,1),x(:,2))
- EnergyStoredInCompressedGas(x(:,1),x(:,2)) - (1./2).*m_s.*(max(vel_s)).^2 )
./(m_p + 4.*I_p.*pi.^2.*R.^2)));
figure
plot(x(:,2),vel_p)
xlabel('x_p')
ylabel('v_p')
title('Velocity of Pellet in Barrel as a Function of Pellet Position in Barrel')
shg

pause

%Next make a plot for (x_p,P(x_p))
Pres = ( (P_o.*V_o)./(V_o - pi.*r_piston.^2.*x(:,1) + (pi./4).*c.^2.*x(:,2)) ) - P_o;
figure
plot(x(:,2),Pres)
xlabel('x_p')
ylabel('P(x_p)')
title('Internal Pressure Within Spring-Air Pellet Gun as a Function of Pellet Position in Barrel')
shg

%Find and display the muzzle velocity ( m/s )
[p8,p8s,p8mu] = polyfit(x(:,2),vel_p,8);
muzzle_velocity = polyval(p8,barrel_length,p8s,p8mu)

%Find and display the muzzle energy ( J )
muzzle_energy = (m_p + 4.*I_p.*pi.^2.*R.^2).*(muzzle_velocity).^2

%Find and display the Maximum Internal Pressure Within Spring-Air Pellet Gun ( Pa )
max_pressure = max(abs(Pres))
%sub-function Library for this file

%Preliminary Functions

function f1 = SpringEnergy(x_s)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f1 = k.*(x_s_f.*x_s - (1./2).*x_s.^2);

function f2 = SimpleGravitationalFriction(x_p)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f2 = m_p.*g.*mu.*x_p;

function f4 = LognitudinalRiflingFriction(x_s,x_p)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f4 = (( (4.*I_p.*R.*pi.*mu.*P_o.*V_o)./(m_p.*c.*(cos(R_angle) + mu.*sin(R_angle))) )
.*log(V_o - pi.*r_piston.^2.*x_s + (1./4).*pi.*c.^2.*x_p) - ( (I_p.*R.*pi.^2.*mu.*c.*P_o.*x_p)
./(m_p.*(cos(R_angle) + mu*sin(R_angle))) ));

function f5 = ElasticCompressionFriction(x_p)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f5 = f_e.*mu.*x_p;

function f6 = BoreFriction(x_s,x_p)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f6 = SimpleGravitationalFriction(x_p) + ElasticCompressionFriction(x_p) + LognitudinalRiflingFriction(x_s,x_p);

function f7 = EnergyStoredInCompressedGas(x_s,x_p)
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
f7 = P_o.*(( (V_o.^2)./(V_o - pi.*r_piston.^2.*x_s + (pi./4).*c.^2.*x_p) ) - V_o - pi.*r_piston.^2.*x_s
+ (pi./4).*c.^2.*x_p);

%Primary ODE system of two linear equations function
function velocity = problem_func(t,x )
global m_p m_s r_piston c R f_e mu I_p g P_o V_o k x_s_f x_s_o R_angle;
%Initialize for a system of ODEs
velocity = zeros(2,1);
%system of two first order ODEs
velocity(1)=(x(1)<(x_s_f-x_s_o)).*(k.*(x_s_f - (x(1) + ((x_s_f-x_s_o)-x(1)).*(x(1) >= (x_s_f-x_s_o))))
- ( (pi.*r_piston^2.*P_o.*V_o)./(V_o - pi.*r_piston.^2.*(x(1) + ((x_s_f-x_s_o)-x(1)).*(x(1)
>= (x_s_f-x_s_o))) + (pi./4).*c.^2.*x(2)) ) + pi.*r_piston.^2.*P_o).*((t)./(m_s));
velocity(2)=sqrt( (SpringEnergy((x(1) + ((x_s_f-x_s_o)-x(1)).*(x(1) >= (x_s_f-x_s_o)))) - BoreFriction((x(1)
+ ((x_s_f-x_s_o)-x(1)).*(x(1) >= (x_s_f-x_s_o))),x(2)) - EnergyStoredInCompressedGas((x(1)
+ ((x_s_f-x_s_o)-x(1)).*(x(1) >= (x_s_f-x_s_o))),x(2)) - (1./2).*m_s.*(max(velocity(1))).^2 )
./(m_p + 4.*I_p.*pi.^2.*R.^2) );
