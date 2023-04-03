%% MEMBERSHIP FUNCTION 'POSITION' PLOT %%

%%% RUN AS A SECTION TO GET THE GRAPH %%%
clc
clear


hold on
figure(1)
% large n
x=-20:-15;
y=1;
hold on
plot(x,y*ones(size(x)),'y')
hold off

hold on
x=-15:-9;
y=-1/6*x-3/2;
plot(x,y,'y')
hold off

hold on
% negative
x=-12:-6;
y=1;
plot(x,y*ones(size(x)),'g')
hold off

hold on
x=-18:-12;
y=1/6*x+3;
plot(x,y,'g')
hold off
hold on
x=-6:0;
y=-1/6*x;
plot(x,y,'g')
hold off

hold on
% zero
x=-3:3;
y=1;
plot(x,y*ones(size(x)),'r')
hold off

hold on
x=-9:-3;
y=1/6*x+3/2;
plot(x,y,'r')
hold off
hold on
x=3:9;
y=-1/6*x+3/2;
plot(x,y,'r')
hold off

hold on
% positive
x=6:12;
y=1;
plot(x,y*ones(size(x)),'b')
hold off

hold on
x=0:6;
y=1/6*x;
plot(x,y,'b')
hold off
hold on
x=12:18;
y=-1/6*x+3;
plot(x,y,'b')
hold off

hold on
% large positive
x=15:20;
y=1;
plot(x,y*ones(size(x)),'m')
hold off

hold on
x=9:15;
y=1/6*x-3/2;
plot(x,y,'m')
hold off


title('Position Membership Functions (MF)')
xlabel('Position (m)')
ylabel('Membership Index')
ylim([0 1.2]);
grid on
hold off

%% MEMBERSHIP FUNCTION 'VELOCITY' PLOT %%

clear;
close all

hold on
figure(2)
% large n
x=-10:-7;
y=1;
plot(x,y*ones(size(x)),'y')
hold on

x=-7:-5;
y=-1/2*x-5/2;
plot(x,y,'y')
hold on

% negative
x=-7:-5;
y=1/2*x+7/2;
plot(x,y,'g')
hold on

x=-5:-3;
y=1;
plot(x,y*ones(size(x)),'g')
hold on

x=-3:-1;
y=-1/2*x-1/2;
plot(x,y,'g')

% zero
x=-3:-1;
y=1/2*x+3/2;
plot(x,y,'r')
hold on

x=-1:1;
y=1;
plot(x,y*ones(size(x)),'r')
hold on

x=1:3;
y=-1/2*x+3/2;
plot(x,y,'r')

% positive
x=1:3;
y=1/2*x-1/2;
plot(x,y,'b')
hold on

x=3:5;
y=1;
plot(x,y*ones(size(x)),'b')
hold on

x=5:7;
y=-1/2*x+7/2;
plot(x,y,'b')

% large positive
x=5:7;
y=1/2*x-5/2;
plot(x,y,'m')

x=7:10;
y=1;
plot(x,y*ones(size(x)),'m')
hold on

title('Velocity Error Membership Functions (MF)')
xlabel('Velocity (m/s)')
ylabel('Membership Index')
ylim([0 1.2]);
grid on
hold off

%% MEMBERSHIP FUNCTION 'BOAT THRUST' PLOT %%

hold on
figure(3)
% large n
x=-300:-250;
y=1;
plot(x,y*ones(size(x)),'y')
hold on

x=-250:-150;
y=-1/100*x-1.5;
plot(x,y,'y')
hold on

% negative
x=-200:-100;
y=1;
plot(x,y*ones(size(x)),'g')
hold on

x=-300:-200;
y=1/100*x+3;
plot(x,y,'g')
hold on
x=-100:0;
y=-1/100*x;
plot(x,y,'g')

% zero
x=-50:50;
y=1;
plot(x,y*ones(size(x)),'r')
hold on

x=-150:-50;
y=1/100*x+1.5;
plot(x,y,'r')
hold on
x=50:150;
y=-1/100*x+1.5;
plot(x,y,'r')
 
%positive
x=100:200;
y=1;
plot(x,y*ones(size(x)),'b')
hold on

x=0:100;
y=1/100*x;
plot(x,y,'b')
hold on

x=200:300;
y=-1/100*x+3;
plot(x,y,'b')

% large positive
x=250:300;
y=1;
plot(x,y*ones(size(x)),'m')
hold on

x=150:250;
y=1/100*x-1.5;
plot(x,y,'m')


title('Thrust Membership Functions (MF)')
xlabel('Thrust (N)')
ylabel('Membership Index')
ylim([0 1.2]);
grid on

hold off

%% INFERENCE RULE MATRIX IRM %%

% Linguistic Values for Position and Velocity Error

dx = 5;

mu_x = membership_function_x(dx);
LN = mu_x(1,1);
N = mu_x(1,2);
Z = mu_x(1,3);
P = mu_x(1,4);
LP = mu_x(1,5);

mu_vel = membership_function_vel(dx);
ln = mu_vel(1,1);
n = mu_vel(1,2);
z = mu_vel(1,3);
p = mu_vel(1,4);
lp = mu_vel(1,5);


X = [LN;N;Z;P;LP]; % Position %
V = [ln n z p lp]; % Velocity Error %

irm = X*V % Inferenece Rule Matrix %
A = irm(4,3) % Recall specific matrix element (Row, Col) %
disp (irm)
t = 0:0.1:1; %s



%% MEMBERSHIP FUNCTION 'VELOCITY' %%

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

%% MEMBERSHIP FUNCTION 'POSITION' %%

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


%% MEMBERSHIP FUNCTION 'THRUST' %%


function mu_thrust = membership_function_thrust(dx)
%Membership functions for the linguistic values of output thrust
%mu_thrust contains the following values (in order):
%  mu_thrust = [LN,N,Z,P,LP]

% large negative (LN)
if dx <= -250
    mu_thrust(1) = 1;
elseif dx > -250 && dx <= -150
    mu_thrust(1) = -1/100*dx-1.5;
else
    mu_thrust(1) = 0;
end

% negative (N)
if dx >= -200 && dx <= -100
    mu_thrust(2) = 1;
elseif dx >= -300 && dx < -200
    mu_thrust(2) = 1/100*dx+3;
elseif dx > -100 && dx <= 0
    mu_thrust(2) = -1/100*dx;
else
    mu_thrust(2) = 0;
end

% zero
if dx >= -50 && dx <= 50
    mu_thrust(3) = 1;
elseif dx >= -150 && dx < -50
    mu_thrust(3) = 1/100*dx+1.5;
elseif dx > 50 && dx <= 150
    mu_thrust(3) = -1/100*dx+1.5;
else
    mu_thrust(3) = 0;
end

% positive
if dx >= 100 && dx <= 200
    mu_thrust(4) = 1;
elseif dx >= 0 && dx < 100
    mu_thrust(4) = 1/100*dx;
elseif dx > 200 && dx <= 300
    mu_thrust(4) = -1/100*dx+3;
else
    mu_thrust(4) = 0;
end

% large positive
if dx >= 250
    mu_thrust(5) = 1;
elseif dx > 150 && dx <= 250
    mu_thrust(5) = 1/100*dx-1.5;
else
    mu_thrust(5) = 0;
end

end


%% DEFUZZIFICATION %%

function T = defuzz(irm_fuzz)

%trapezoid base sizes
LB = 300;
SB = 100;

%area calulcations
area_mat = zeros(5,5);

a = 1;
b = 1;

while a < 6
    while b < 6
        area_mat(a,b) = LB*irm_fuzz(a,b) - (((irm_fuzz(a,b))^2)/2)*(LB - SB);
        b = b + 1;
        
    end
    b = 1;
    a = a + 1;
end

%centroids
CLN = -300;
CN = -150;
CZ = 0;
CP = 150;
CLP = 300;

cent_mat = [CLP,CLP,CLP,CP, CZ;...
            CLP,CP, CP, CZ, CN;...
            CLP,CP, CZ, CN, CLN;...
            CP, CZ, CN, CN, CLN;...
            CZ, CN, CLN,CLN,CLN];

%summation
%centroid*area
numerator = 0;
a = 1;
b = 1;
while a < 6
    while b < 6
        numerator = numerator + area_mat(a,b)*cent_mat(a,b);
        b = b + 1;
        
    end
    b = 1;
    a = a + 1;
end

%denominator
a = 1;
b = 1;
denom = 0;
while a < 6
    while b < 6
        denom = denom + area_mat(a,b);
        
        b = b + 1;
        
    end
    b = 1;
    a = a + 1;
end

T = numerator/denom;

end
