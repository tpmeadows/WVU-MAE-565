%Route-finding GA

clear
close all
clc

% bike.cap = 15; %# of packages
% bike.v = 45; %mph
% bike.mpg = 200; %mpg
% 
% car.cap = 35;
% car.v = 40;
% car.mpg = 24;
% 
% truck.cap = 50;
% truck.v = 35;
% truck.mpg = 9;

v = 45; %mph
cap = 20; %carrying capacity

n = 100; %pop size
moves = [0 1;1 0;0 -1;-1 0];
m = 0.1*n; %# of mutated individuals
c = 0.5*n; %# of reproducing individuals
L0 = 8*v+2; %max # of stops
pop(1:n) = struct('path',zeros(L0,2),'turns',zeros(L0,1),'distance',0,'delivered',0,'checklist',0);

%Criteria
chk0 = 150; %# of visited houses
d0 = 100; %distance traveled
deliver0 = 150; %# of deliveries

%%%%%%%%%       Initialization       %%%%%%%%%
EF = zeros(n,1);
for i = 1:n
    %startveh = randi(3);
    chk = 0; %counter for first visits
    tr = 8; %max time, hours
    delivered = -1; %packages delivered. starting @ -1 to remove origin from route
    stops = [0 0];
    %a = []; %storage variable
    direction = zeros(L0,1);
    total_distance = 0;
    j = 1;
    while tr >= 0
       %for directions, N=1, E=2, S=3, W=4
       %k vector represents the set of legal moves
%        if chk == 0
%            chk = chk0;
%            a = [a;stops];
%            stops = stops(end,:);
%        end

       if stops(j,1) == 0
           if stops(j,2) == 0
               k = [1 2];
               carry = cap;
           elseif stops(j,2) == 9
               k = [2 3];
           else
               k = [1 2 3];
           end
       elseif stops(j,1) == 9
           if stops(j,2) == 0
               k = [1 4];
           elseif stops(j,2) == 9
               k = [3 4];
           else
               k = [1 3 4];
           end
       elseif stops(j,2) == 0
           k = [1 2 4];
       elseif stops(j,2) == 9
           k = [2 3 4];
       else
           k = [1 2 3 4];
       end

       l = length(k);
       direction(j) = k(randi(l));
       new_stop = stops(j,:) + moves(direction(j),:);
       if isempty(find(stops(:,1) == new_stop(1) & stops(:,2) == new_stop(2),1)) == true ...
               && carry ~= 0
           chk = chk + 1;
           carry = carry - 1;
           delivered = delivered + 1;
       end
       stops(j+1,:) = new_stop;

       j = j+1;
       d = 1/v*sqrt((stops(j,1)-stops(j-1,1))^2 + (stops(j,2)-stops(j-1,2))^2);
       total_distance = total_distance + d;
       tr = tr - d;
    end
    pop(i).path = stops;%[a;stops];
    pop(i).turns = direction;
    pop(i).distance = total_distance;
    pop(i).delivered = delivered;
    pop(i).checklist = chk;
    EF(i,1) = 0.5*(chk/chk0) + 0.3*(total_distance/d0) + 0.2*(delivered/deliver0);
end

%%%%%%%%% Selection (Roulette Wheel) %%%%%%%%%
TF = sum(EF);
q = zeros(n,1);
for i = 1:n
    q(i) = sum(EF(1:i))/TF;
end
pop0(1:n) = struct('path',zeros(L0,2),'turns',zeros(l,1),'distance',0,'delivered',0,'checklist',0);
I = ones(n,1);
for i = 1:n
    r = rand;
    a = abs(q - r*I);
    j = find(a == min(a));
    pop0(i) = pop(j);
end

%%%%%%%%%     Genetic Operators     %%%%%%%%%
ind_m = randi(n,m,1); %indices of mutating individuals
pop_m = pop0(ind_m);

for i = 1:m
    r = randi([round(0.75*n) n]);
    a = pop_m.path(r,:);
    if a(i,1) == 0
        if a(i,2) == 0
            k = [1 2];
            carry = cap;
        elseif a(i,2) == 9
            k = [2 3];
        else
            k = [1 2 3];
        end
    elseif a(i,1) == 9
        if a(i,2) == 0
            k = [1 4];
        elseif stops(j,2) == 9
            k = [3 4];
        else
            k = [1 3 4];
        end
    elseif a(i,2) == 0
        k = [1 2 4];
    elseif a(i,2) == 9
        k = [2 3 4];
    else
        k = [1 2 3 4];
    end
    %%%%% WIP %%%%%
end
