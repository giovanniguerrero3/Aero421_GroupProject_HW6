function [ c21_matrix ] = rot321( euler )
% Rotation is from ECI to LVLH
%*** Keep all angles in RADIANS ***

Cx = [1 0 0;
    0 cos(euler(1)) sin(euler(1));
    0 -sin(euler(1)) cos(euler(1))];

Cy = [cos(euler(2)) 0 -sin(euler(2));
    0 1 0;
    sin(euler(2)) 0 cos(euler(2))];

Cz = [cos(euler(3)) sin(euler(3)) 0;
    -sin(euler(3)) cos(euler(3)) 0;
    0 0 1];

c21_matrix = Cx*Cy*Cz;

end

