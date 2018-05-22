function [R,V] = coes2rv(h,mu,ecc,theta,RAAN,inc,w)
%COES2RV
%All angles inputted as degrees

RAAN = RAAN*pi/180;
i = inc*pi/180;
w = w*pi/180;

Rpqw = h^2/(mu*(1+ecc*cosd(theta))).*[cosd(theta);sind(theta);0];
Vpqw = mu/h.*[-sind(theta); ecc+cosd(theta); 0];

Q = [-sin(RAAN)*cos(i)*sin(w)+cos(RAAN)*cos(w), cos(RAAN)*cos(i)*sin(w)+sin(RAAN)*cos(w), sin(i)*sin(w);
    -sin(RAAN)*cos(i)*cos(w)-cos(RAAN)*sin(w), cos(RAAN)*cos(i)*cos(w)-sin(RAAN)*sin(w), sin(i)*cos(w);
    sin(RAAN)*sin(i), -cos(RAAN)*sin(i), cos(i)];

%Must use the inverse of Q to go from perifocal to ECI

R = Q'*Rpqw;
V = Q'*Vpqw;

end