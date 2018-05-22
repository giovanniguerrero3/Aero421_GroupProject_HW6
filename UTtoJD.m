function [ J0,JD ] = jd_calc( y,m,d,UT )
% Input: year [****], month [**], day [**], time (UT)
% Output: J0, JD

J0=(367*y) - floor(7*(y+floor((m+9)/12))/4) + floor((275*m)/9) + d + 1721013.5;
JD=J0+(UT/24);

end

