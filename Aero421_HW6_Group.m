% Giovanni Guerrero, Roberto Azarte, Chris Barta
% Aero 421 Group Assignment 6

clear; close all; clc
format long g
%% Set Up

%COES
global mu I com 
mu =398600;

period = 100*60; %period in seconds (100 mins);

h = 53335.2; %km^2/s (angular momentum)
ecc = 0; %eccentricity)
RAAN = 0; %deg
inc = 98.43; %deg (inclination)
w = 0; %deg (omega (argument of periapsis))
ta = 0; %deg (true anomaly)

[Ri,Vi] = coes2rv(h,mu,ecc,ta,RAAN,inc,w);

w_bi = [0; -0.001047; 0]; % (rad/s)initial in body frame

rho = 1.139; %kg/m^3
a_panel = 3*2; %[m^2]
a_body = 2^2; %[m^2]
Cd = 2.1; %coefficient of drag

eps = [0 0 0];
eta = 1;
    
q = [eps eta]'; %"not initial quaternion"

%simplified density, drag, solar pressure, and earth magnetic field models

% Spacecraft properties 
c_m = 500; %kg
c_sx = 2; %m, side of cube
c_sy = 2;
c_sz = 2;

c_dz = 0;
c_dy = 0.234375;
c_dx = 0.234375;

%sensors
s_m = 100; %(kg) mass of sensor
s_sz = 1; %length of sensor
s_sx = 0.25; %length in x
s_sy = 0.25; %length in y

s_dx = 1.265625; % center of sensor to center of bus
s_dy = 1.265625;
s_dz = 0; 

% Solar Panels (accounting for both)
p_x = 2; %x and y solar panel face
p_y = 3;
p_z = 0.05;
p_m = 20; %20kg each

p_dx = 2.5; % center of panels to center of bus
p_dy = 0.23437;
p_dz = 2.5;

%% Problem 1 Properties of Mass and Inertia
disp('----------------------------')
disp('Problem 1')
disp('----------------------------')

%total mass of s/c
m = c_m+s_m+2*p_m;


%Calculating inertias relative to total COM
%Bus Cube
c_Ix = (c_m/12)*(c_sy^2+c_sz^2);
c_Iy = (c_m/12)*(c_sx^2+c_sz^2);
c_Iz = (c_m/12)*(c_sy^2+c_sx^2);
I_c = [c_Ix+c_m*c_dx^2 0 0; 0 c_Iy+c_m*c_dy^2 0; 0 0 c_Iz+c_m*c_dz^2];

%Sensors
s_Ix = (s_m/12)*(s_sy^2+s_sz^2);
s_Iy = (s_m/12)*(s_sx^2+s_sz^2);
s_Iz = (s_m/12)*(s_sy^2+s_sx^2);
I_s = [s_Ix+s_m*s_dx^2 0 0;0 s_Iy+s_m*s_dy^2 0;0 0 s_Iz+s_m*s_dz^2];

%Solar panels
p_Ix = (p_m/12)*(p_y^2+p_z^2);
p_Iy = (p_m/12)*(p_z^2+p_x^2);
p_Iz = (p_m/12)*(p_y^2+p_x^2);
I_p = 2*[p_Ix+p_m*p_dx^2 0 0;0 p_Iy+p_m*p_dy^2 0;0 0 p_Iz+p_m*p_dz^2];


%total Inertia of s/c

I = I_c+I_p+I_s;

disp(['Mass of spacecraft is ',num2str(m),' kg'])
disp(' ')
disp(['Inertia of spacecraft is ',mat2str(I),' kg*m^2']);
disp(' ')
disp('----------------------------')

%Note: Origin of body frame is at center of bus

%Center of mass (relative to center of bus(done on board))
 com = [0,0,0.23437];
 com = [0 0 0; 0 0 0; 0 0 0.23437]; %is this the 3D one?
 
%Geometric Center (relative to center of bus(done on paper))


%% Problem 2

disp('Problem 2')
disp('----------------------------')
% 
% %part a)
% init_theta = atan(Vi(2)/Vi(3)); %[rad] angle at the beginning of orbit between Xb and Zeci
% 
% 
% %initial angular momentum
% h_i=I*w_bi; %[kg*m^2/s] (3x1)
% 
% %initial angular velocity (cross)
% w_b_cross=cross_matrix(w_bi); %rad/s
% 
% %kinematics
% vect = [ init_theta 0 0];
% [c21_matrix] = rot321(vect);
% wbi_x = cross_matrix(w_bi);
% 
% Cbi_dot = -wbi_x*c21_matrix;
% 
%part b)
Fd_panel = (1/2)*rho*Cd*a_panel*norm(Vi)^2;
Fd_body = (1/2)*rho*Cd*a_body*norm(Vi)^2;

% %part c)
% 
% %part d)
% % r is the r vector from the center of the earth to the center of the body?
% Tg = (3*mu)/(norm(r)^5)*cross(r,I*r); % Gravity Gradient Torque, taken from notes on April 10
% 
% %part e)


%% Propogate and Plot

% Plotting 1 period and marking starting point
statei = [Ri Vi];
tspan = [0 period];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

[t_orbit,states_orbit] = ode45(@propogate_func, tspan, statei, options, mu);

rstart = [states_orbit(end,1),states_orbit(end,2),states_orbit(end,3)]; %basically Ri

% R and V vectors changing with time throughout one period. In ECI Frame.
R_state = [states_orbit(:,1),states_orbit(:,2),states_orbit(:,3)];
V_state = [states_orbit(:,4),states_orbit(:,5),states_orbit(:,6)];
time_state = t_orbit;

figure 
hold on 
plot3(rstart(1),rstart(2),rstart(3),'x','Color','r','linewidth',5)
plot3(states_orbit(:,1), states_orbit(:,2), states_orbit(:,3), 'Color', 'r')
grid on
view(3)



%Adding Earth to Plot
grs80 = referenceEllipsoid('grs80','km');
ax = axesm('globe','Geoid',grs80,'Grid','off');
ax.Position = [0 0 1 1];
axis equal off
view(3)

load topo
geoshow(topo,topolegend,'DisplayType','texturemap')
demcmap(topo)
land = shaperead('landareas','UseGeoCoords',true);
plotm([land.Lat],[land.Lon],'Color','black')
rivers = shaperead('worldrivers','UseGeoCoords',true);
plotm([rivers.Lat],[rivers.Lon],'Color','blue')


%% March 20, 2018 at 12 UTC: Julian Date and Epoch Location  
% Calculate Julian Date of epoch
% Epoch = 17325.34312598 
y=2017;              % year
m=11;                % month
d=21;                % day
UT=.34312598*24;     % UT
[ J01, JD1 ] = UTtoJD( y,m,d,UT );

% Mission Start Date: Spring Equinox: March 20, 2018, UT=12 (noon)
[ J0start, JDstart ] = UTtoJD( 2018,3,20,12 );

dJD1=JDstart-JD1; 
JD2sec=.864/0.00001;    % conversion factor from julian days to seconds
dt1=dJD1*JD2sec;        % time from epoch to start date (March 20, 2018)

% Locating Epoch (start)
%COES have to be in radians for Locate function. 
h = 53335.2; %km^2/s (angular momentum)
ecc = deg2rad(ecc); %rad eccentricity)
RAAN = deg2rad(RAAN); %rad
inc = deg2rad(inc); %rad (inclination)
w = deg2rad(w); %rad (omega (argument of periapsis))
ta = deg2rad(ta); %rad (true anomaly)
% Object Location at Epoch
[ ~, ~, ri_epoch, vi_epoch ] = Locate( ecc, inc, RAAN, w, ta, h );

%% Functions

%function to be integrated using ode 45
function state_out=ode_funct(t,state,T,I)

%euler angles
eu_ang=[state(4);state(5);state(6)];

%rotation matrix
C21=rotation_mtrx(eu_ang);

%quaternion
q=quat(C21);

%angular velocity
ang_vel=[state(1);state(2);state(3)];

%angular velocity rates
ang_vel_rates=euler_eqns(ang_vel,I,T);

%euler angle rates
eu_ang_rates=euler_rates(eu_ang,ang_vel);

%quaternion rates
q_rates=quat_rates(q,ang_vel);

%output state vector
state_out=[ang_vel_rates;eu_ang_rates;q_rates];

end

%calculate rotation matrix
function C21=rotation_mtrx(ang)

Cx=[1 0 0; 0 cos(ang(1,1)) sin(ang(1,1)); 0 -sin(ang(1,1)) cos(ang(1,1))];
Cy=[cos(ang(2,1)) 0 -sin(ang(2,1)); 0 1 0; sin(ang(2,1)) 0 cos(ang(2,1))];
Cz=[cos(ang(3,1)) sin(ang(3,1)) 0; -sin(ang(3,1)) cos(ang(3,1)) 0; 0 0 1];

%3-2-1 rotation matrix
C21=Cx*Cy*Cz;

end

%calculate quaternion
function q=quat(C21)

%eta
eta=sqrt(trace(C21)+1)/2;

%epsilon components
e1=(C21(2,3)-C21(3,2))/(4*eta);
e2=(C21(3,1)-C21(1,3))/(4*eta);
e3=(C21(1,2)-C21(2,1))/(4*eta);

%epsilon
eps=[e1;e2;e3];

%quaternion
q=[eps;eta];

end

%angular velocity rates
function ang_vel_rates=euler_eqns(w,I,T)

%euler equations
w_dx=(T(1,1)+w(2,1)*w(3,1)*(I(2,2)-I(3,3)))/I(1,1);
w_dy=(T(2,1)+w(1,1)*w(3,1)*(I(3,3)-I(1,1)))/I(2,2);
w_dz=(T(3,1)+w(1,1)*w(2,1)*(I(1,1)-I(2,2)))/I(3,3);

%equations as column vector
ang_vel_rates=[w_dx;w_dy;w_dz];

end

%calculates euler angle rates
function eu_ang_rates=euler_rates(angles,w)

eu_ang_rates=...
    [1 sin(angles(1))*tan(angles(2)) cos(angles(1))*tan(angles(2)); 
    0 cos(angles(1)) -sin(angles(1)); 
    0 sin(angles(1))/cos(angles(2)) cos(angles(1))/cos(angles(2))]*w;

end

%calculates eta dot & epsilon dot
function q_rates=quat_rates(q,w)

%epsilon cross
eps_x=[0 -q(3,1) q(2); q(3,1) 0 -q(1,1); -q(2,1) q(1,1) 0];

%eta dot
eta_d=-0.5*q(1:3,1)'*w;

%epsilon dot
eps_d=0.5*(q(4,1)*eye(3)+eps_x)*w;

%quaternion rate function
q_rates=[eps_d;eta_d];

end

%solar pressure torque
function T_s=sol_press_torq(s_vect)

%top surfaces (normal vectors)
n_top=[0; 0; -1];

%bottom surfaces (normal vectors)
n_bot=[0; 0; 1];

%front and back
n_f=[1; 0; 0];
n_b=[-1; 0; 0];

%sides
n_l=[0; -1; 0];
n_r=[0; 1; 0];

%normal vector, vector
n=[n_top n_top n_top n_bot n_bot n_bot n_f n_b n_l n_r];

%position vectors
a=0.234375; %m
b=0.7656; %m
c=1.7656; %m
d=2.5; %m
top_c=[0;0;-c];
bot_c=[0;0;b];
front_c=[1;0;-a];
back_c=[-1;0;-a];
l_side_c=[0;-1;-a];
r_side_c=[0;1;-a];
l_sol_p_top=[0;-d;-a];
l_sol_p_bot=[0;-d;-a];
r_sol_p_top=[0;d;-a];
r_sol_p_bot=[0;d;-a];

%position relative to com, vector
rho=[top_c l_sol_p_top r_sol_p_top bot_c l_sol_p_bot r_sol_p_bot...
    front_c back_c l_side_c r_sidec];

%photon momentum
p=4.5e-6; %N*m^-2

%spacecraft areas
A_c=4; %m^2
A_p=6; %m^2

A=[A_c A_p A_p A_c A_p A_p A_c A_c A_c A_c];

nds=dot(n,s);

if nds>=0

T_s=cross(rho,-p.*A*nds*s);

else
    T_s=0;

end
