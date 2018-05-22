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

w_bi = [0 -0.001047 0]; % (rad/s)initial in body frame

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

%part a)

%part b)

%part c)

%part d)
% r is the r vector from the center of the earth to the center of the body?
Tg = (3*mu)/(norm(r)^5)*cross(r,I*r); % Gravity Gradient Torque, taken from notes on April 10

%part e)


%% Propogate and Plot

% Plotting 1 period and marking starting point
statei = [Ri Vi];
tspan = [0 period];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

[t,statenew] = ode45(@propogate_func, tspan, statei, options, mu);

rstart = [statenew(end,1),statenew(end,2),statenew(end,3)]; %basically Ri

figure 
hold on 
plot3(rstart(1),rstart(2),rstart(3),'x','Color','r','linewidth',5)
plot3(statenew(:,1), statenew(:,2), statenew(:,3), 'Color', 'r')
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













