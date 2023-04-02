clear;
%crisp input
vel_input = -6;


if vel_input < -7
    fuzzy_vel = [1,0,0,0,0];
elseif vel_input < -5
    LN = (-1/2*vel_input)-5/2;
    N = (1/2*vel_input)+7/2;
    fuzzy_vel = [LN,N,0,0,0];
elseif vel_input < -3
    fuzzy_vel = [0,1,0,0,0];
elseif vel_input < -1
    N = (-1/2*vel_input)-1/2;
    Z = (1/2*vel_input)+3/2;
    fuzzy_vel = [0,N,Z,0,0];
elseif vel_input < 1
    fuzzy_vel = [0,0,1,0,0];
elseif vel_input < 3
    Z = (-1/2*vel_input)+3/2;
    P = (1/2*vel_input)-1/2;
    fuzzy_vel = [0,0,Z,P,0];
elseif vel_input < 5
    fuzzy_vel = [0,0,0,1,0];
elseif vel_input < 7
    P = (-1/2*vel_input)+7/2;
    L = (1/2*vel_input)-5/2;
    fuzzy_vel = [0,0,0,P,L];
else
    fuzzy_vel = [0,0,0,0,1];
end