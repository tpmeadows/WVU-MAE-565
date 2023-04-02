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