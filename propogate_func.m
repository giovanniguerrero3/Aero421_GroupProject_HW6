function [ dyydt ] = propogate_func( tspan,state, mu )

r = norm([state(1) state(2) state(3)]);
dxdt = -mu*state(1)/(r^3);
dydt = -mu*state(2)/(r^3);
dzdt = -mu*state(3)/(r^3);

dyydt = [state(4); state(5); state(6); dxdt; dydt; dzdt];


end

