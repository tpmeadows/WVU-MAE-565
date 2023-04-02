%% Functions

function mu_x = membership_function_x(dx)
%Membership functions for the linguistic values of position error
%mu_x contains the following values (in order):
%  mu_x = [LN,N,Z,P,LP]

% large negative (LN)
if dx <= -15
    mu_x(1) = 1;
elseif dx > -15 && dx <= -9
    mu_x(1) = -1/6*dx-3/2;
else
    mu_x(1) = 0;
end

% negative (N)
if dx >= -12 && dx <= -6
    mu_x(2) = 1;
elseif dx >= -18 && dx < -12
    mu_x(2) = 1/6*dx+3;
elseif dx > -6 && dx <= 0
    mu_x(2) = -1/6*dx;
else
    mu_x(2) = 0;
end

% zero
if dx >= -3 && dx <= 3
    mu_x(3) = 1;
elseif dx >= -9 && dx < -3
    mu_x(3) = 1/6*dx+3/2;
elseif dx > 3 && dx <= 9
    mu_x(3) = -1/6*dx+3/2;
else
    mu_x(3) = 0;
end

% positive
if dx >= 6 && dx <= 12
    mu_x(4) = 1;
elseif dx >= 0 && dx < 6
    mu_x(4) = 1/6*dx;
elseif dx > 12 && dx <= 18
    mu_x(4) = -1/6*dx+3;
else
    mu_x(4) = 0;
end

% large positive
if dx >= 15
    mu_x(5) = 1;
elseif dx > 9 && dx <= 15
    mu_x(5) = 1/6*dx-3/2;
else
    mu_x(5) = 0;
end

end

function mu_vel = membership_function_vel(vel_input)

% large negative (LN)
if vel_input < -7
    mu_vel = [1,0,0,0,0];

elseif vel_input < -5
    LN = (-1/2*vel_input)-5/2;
    N = (1/2*vel_input)+7/2;
    mu_vel = [LN,N,0,0,0];
    
% negative (N)
elseif vel_input < -3
    mu_vel = [0,1,0,0,0];
    
elseif vel_input < -1
    N = (-1/2*vel_input)-1/2;
    Z = (1/2*vel_input)+3/2;
    mu_vel = [0,N,Z,0,0];

% zero (Z)
elseif vel_input < 1
    mu_vel = [0,0,1,0,0];
    
elseif vel_input < 3
    Z = (-1/2*vel_input)+3/2;
    P = (1/2*vel_input)-1/2;
    mu_vel = [0,0,Z,P,0];
    
%positive (P)   
elseif vel_input < 5
    mu_vel = [0,0,0,1,0];
    
elseif vel_input < 7
    P = (-1/2*vel_input)+7/2;
    L = (1/2*vel_input)-5/2;
    mu_vel = [0,0,0,P,L];
    
%large positive (LP)
else
    mu_vel = [0,0,0,0,1];
    
end

end

function mu_thrust = membership_function_thrust(dx)
%Membership functions for the linguistic values of output thrust
%mu_thrust contains the following values (in order):
%  mu_thrust = [LN,N,Z,P,LP]

% large negative (LN)
if dx <= -150000
    mu_thrust(1) = 1;
elseif dx > -150000 && dx <= -90000
    mu_thrust(1) = -1/60000*dx-3/2;
else
    mu_thrust(1) = 0;
end

% negative (N)
if dx >= -120000 && dx <= -60000
    mu_thrust(2) = 1;
elseif dx >= -180000 && dx < -120000
    mu_thrust(2) = 1/60000*dx+3;
elseif dx > -60000 && dx <= 0
    mu_thrust(2) = -1/60000*dx;
else
    mu_thrust(2) = 0;
end

% zero
if dx >= -30000 && dx <= 30000
    mu_thrust(3) = 1;
elseif dx >= -90000 && dx < -30000
    mu_thrust(3) = 1/60000*dx+3/2;
elseif dx > 30000 && dx <= 90000
    mu_thrust(3) = -1/60000*dx+3/2;
else
    mu_thrust(3) = 0;
end

% positive
if dx >= 60000 && dx <= 120000
    mu_thrust(4) = 1;
elseif dx >= 0 && dx < 60000
    mu_thrust(4) = 1/60000*dx;
elseif dx > 120000 && dx <= 180000
    mu_thrust(4) = -1/60000*dx+3;
else
    mu_thrust(4) = 0;
end

% large positive
if dx >= 150000
    mu_thrust(5) = 1;
elseif dx > 90000 && dx <= 150000
    mu_thrust(5) = 1/60000*dx-3/2;
else
    mu_thrust(5) = 0;
end

end
