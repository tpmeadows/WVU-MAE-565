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