function [ rPER_vec, vPER_vec, rECI_vec, vECI_vec ] = Locate( ecc, inc, RAAN, w, theta, h )
% Input: eccentricity [rad], inclination [rad], RAAN [rad], argument of 
% perigee [rad], true anomaly [rad], altitude at perigee [km]
% Output: perifocal r (p-direction)[km] and v at perigee (q-direction)[km/s] vectors
% at perigee, geocentric r [km] and v [km/s] vectors at perigee


muEarth=398600;      % [km]
rEarth=6378;        % [km]

% 1. Calculate position vector in perifocal coordinates using Eqn (4.45).
rPER_vec=(h^2/muEarth)*(1/(1+(ecc*cos(theta))))*[cos(theta); sin(theta); 0];

% 2. Calculate velocity vector in perifocal coordinates using Eqn (4.46).
vPER_vec=(muEarth/h)*[(-sin(theta)); ecc+cos(theta);0];

% 3. Calculate the matrix QxX of the transformation from perifocal to
% geocentric equatorial coordinates using Eqn (4.49).
QxX=[-sin(RAAN)*cos(inc)*sin(w)+cos(RAAN)*cos(w) , -sin(RAAN)*cos(inc)*cos(w)-cos(RAAN)*sin(w) , sin(RAAN)*sin(inc);...
    (cos(RAAN)*cos(inc)*sin(w))+(sin(RAAN)*cos(w)) , cos(RAAN)*cos(inc)*cos(w)-sin(RAAN)*sin(w) , -cos(RAAN)*sin(inc);...
    sin(inc)*sin(w) , sin(inc)*cos(w) , cos(inc) ];

% 4. Transform rx and vx into the geocentric frame by means of Eqn (4.51).
rECI_vec=QxX*rPER_vec;
vECI_vec=QxX*vPER_vec;

end



