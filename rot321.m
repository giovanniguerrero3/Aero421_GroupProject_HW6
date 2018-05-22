function [ c21_matrix ] = rot321( phi, theta, psy )

% Keep all angles in RADIANS

Cx = [1 0 0;
    0 cos(phi) sin(phi);
    0 -sin(phi) cos(phi)];

Cy = [cos(theta) 0 -sin(theta);
    0 1 0;
    sin(theta) 0 cos(theta)];

Cz = [cos(psy) sin(psy) 0;
    -sin(psy) cos(psy) 0;
    0 0 1];

c21_matrix = Cx*Cy*Cz;

end

