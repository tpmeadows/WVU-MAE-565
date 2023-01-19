
clear
close all
clc

%Optimization Parameters
%route (constant, speed limits), transportation (carry capacity, speed, vary across delivery route, tank volume, mileage)

%Optimization Criteria
%emissions, cost, packages delivered, stops

bike.capacity = [5 20]; %packages
bike.Vol = [3 6]; %gal
bike.mpg = [180 220]; %miles/gal
bike.v = [35 50]; %mph

car.capacity = [20 40]; %packages
car.Vol = [13 16]; %gal
car.mpg = [19 29]; %miles/gal
car.v = [30 45]; %mph

truck.capacity = [40 100]; %packages
truck.Vol = [22 36]; %gal
truck.mpg = [8 10]; %miles/gal
truck.v = [25 40]; %mph

par0 = [1 2 3 4 5];
cri0 = [1 2 3 4 5];
x = (0:10); y = (0:10);
%D = 100; %distance traveled
C = 3.25; % $/gal
em = 8887; %gC02/gal of gasoline
t = 8; %timeframe, hours
l = length(cri0);
cri = cri0;
marg = 0.9*ones(l,1); %percentage agreement

while cri./cri0 >= marg
    
end


function B = GA(A)
    v = A(1);
    Vol = A(2);
    mpg = A(3);
    cap = A(4);
    
    
end

hey zach can you see this???
