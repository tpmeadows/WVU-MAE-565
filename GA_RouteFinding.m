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

n = 200; %pop size
moves = [0 1;1 0;0 -1;-1 0];
m = round(0.10*n); %# of mutating individuals
c = round(0.75*n); %# of reproducing individuals
if rem(c,2) ~= 0
    c = c + 1;
end
L0 = 8*v+2; %max # of stops
pop(1:n) = struct('path',zeros(L0,2),'turns',zeros(L0,1),'distance',0,'delivered',0);

%Criteria
chk0 = 100; %# of visited houses - origin
d0 = 100; %distance traveled
deliver0 = 150; %# of deliveries

figure(1)
title('Delivery Path')

figure(2)
xlabel('Generation')
ylabel('Evaluated Fitness')
ylim([0 1])
grid on
hold on

%%%%%%%%%       Initialization       %%%%%%%%%
EF = zeros(n,1);
for i = 1:n
    %startveh = randi(3);
    tr = 8; %max time, hours
    chk = chk0;
    delivered = 0;
    stops = [0 0];
    a = []; %storage variable
    direction = zeros(L0,1);
    total_distance = 0;
    j = 1;
    while tr >= 0
       %for directions, N=1, E=2, S=3, W=4
       %k vector represents the set of legal moves
       if chk == 0
           chk = chk0;
           a = [a;stops];
           stops = stops(end,:);
       end
       
       x = stops(j,1); y = stops(j,2);
       if x == 0 && y == 0
           carry = cap;
       end
       k = LM(x,y);
       l = length(k);
       direction(j) = k(randi(l));
       new_stop = [x,y] + moves(direction(j),:);
       if isempty(find(stops(:,1) == new_stop(1) & stops(:,2) == new_stop(2),1)) == true ...
               && carry ~= 0 && (x ~= 0 || y ~= 0)
           chk = chk - 1;
           carry = carry - 1;
           delivered = delivered + 1;
       end
       stops(j+1,:) = new_stop;

       j = j+1;
       d = 1/v*sqrt((stops(j,1)-x)^2 + (stops(j,2)-y)^2);
       total_distance = total_distance + d;
       tr = tr - d;
    end
    pop(i).path = stops;%[a;stops];
    pop(i).turns = direction;
    pop(i).distance = total_distance;
    pop(i).delivered = delivered;
    EF(i,1) = 0.2*(total_distance/d0) + 0.8*(delivered/deliver0);
end

M_EF = max(EF);
ind_max = find(EF == M_EF);
ind_max = ind_max(randi(length(ind_max)));
G = 1;
p = 1;
while M_EF <= 0.8

%%%%%%%%% Selection (Roulette Wheel, Elitist) %%%%%%%%%
    
    TF = sum(EF);
    q = zeros(n,1);
    for i = 1:n
        q(i) = sum(EF(1:i))/TF;
    end
    pop0(1:n) = struct('path',zeros(L0,2),'turns',zeros(L0,1),'distance',0,'delivered',0);
    I = ones(n,1);
    for i = 1:n
        r = rand;
        a = abs(q - r*I);
        j = find(a == min(a));
        if i ~= ind_max
            pop0(i) = pop(j);
        else
            pop0(i) = pop(i);
        end
    end

%%%%%%%%%     Genetic Operators     %%%%%%%%%

    %%%%%% Mutation %%%%%%
    ind_m = randi(n,m,1); %indices of mutating individuals
    pop_m = pop0(ind_m);
    %n_m = randi([5 round(L0/8)]); %# of times mutation is applied to an individual
    %for mutations = 1:n_m
        for i = 1:m
            if ind_m ~= ind_max
                r = randi(L0);
                pos_m = pop_m(i).path(r,:);
                turn_m = pop_m(i).turns(r);

                k = LM(pos_m(1),pos_m(2));
                ind = find(k ~= turn_m);
                l = length(ind);
                mv_m = k(ind(randi(l)));
                pop_m(i).turns(r) = mv_m;
                for j = r:L0-1
                    b = pos_m + moves(mv_m,:);
                    if b(1) > 9 || b(1) < 0 || b(2) > 9 || b(2) < 0
                        k = LM(pos_m(1),pos_m(2));
                        ind = find(k ~= turn_m);
                        l = length(ind);
                        mv_m = k(ind(randi(l)));
                        pop_m(i).turns(j) = mv_m;
                        b = pos_m + moves(mv_m,:);
                    end
                    pop_m(i).path(j+1,:) = b;
                    pos_m = b;
                    mv_m = pop_m(i).turns(j+1);
                end
                pop0(ind_m(i)) = pop_m(i);
            end
        end
    %end
    
    %%%%%%%%%    Crossover    %%%%%%%%%
    ind_c = randi(n,c,1); %indices of reproducing individuals
    ind_c(ind_c == ind_max) = [];
    lic = length(ind_c);
    pop_c = pop0(ind_c);
    intvl = [round(L0*0.1) round(L0*0.9)]; %bounds for crossover point

    for i = 1:c/2
        i1 = ind_c(randi(round(lic/2)));
        pos1 = pop0(i1).path;
        turn1 = pop(i1).turns;
        i2 = ind_c(randi([round(lic/2)+1 lic]));
        pos2 = pop0(i2).path;
        turn2 = pop0(i2).turns;

        i_crs = find(~(pos1(:,1)-pos2(:,1)) & ~(pos1(:,2)-pos2(:,2)));
        i_crs = i_crs(i_crs > intvl(1) & i_crs < intvl(2));
        lc = length(i_crs);
        if lc ~= 0
            crosspt = i_crs(randi(lc));
            for j = crosspt:L0
                a = pop0(i1).path(j,:);
                b = pop0(i2).path(j,:);

                pop0(i1).path(j,:) = b;
                pop0(i2).path(j,:) = a;
            end
        end
    end

%%%%%%%%%    Fitness Evaluation     %%%%%%%%%
    for i = 1:n
        chk = chk0;
        carry = cap;
        delivered = 0;
        total_distance = 0;
        stops = [];
        a = 1;
        b = 0;
        for j = 2:L0
            stops = pop0(i).path(a:j,:);
            c_stops = j-b;
            x = stops(c_stops,1); y = stops(c_stops,2);
            if x == 0 && y == 0
                carry = cap;
            end

            if isempty(find(stops(1:c_stops-1,1) == x & stops(1:c_stops-1,2) == y,1)) == true ...
                    && carry ~= 0 && (x ~= 0 || y ~= 0)
                chk = chk - 1;
                carry = carry - 1;
                delivered = delivered + 1;
            end

            if chk == 0
                chk = chk0;
                a = j;
                b = b + j;
            end

            d = 1/v*sqrt((stops(j-1,1)-x)^2 + (stops(j-1,2)-y)^2);
            total_distance = total_distance + d;
        end

        pop0(i).distance = total_distance;
        pop0(i).delivered = delivered;
        EF(i,1) = 0.2*(total_distance/d0) + 0.8*(delivered/deliver0);
    end
    M_EF = max(EF);
    EFavg = mean(EF);
    ind_max = find(EF == M_EF);
    lmax = length(ind_max);
    ind_max = ind_max(randi(lmax));
    pop_f = pop0(ind_max);
    x = pop_f.path(:,1); y = pop_f.path(:,2);
    pop = pop0;

    figure(1)
    plot(x,y)
    xlim([0 9])
    ylim([0 9])
    lgd = strcat('Generation #',string(G));
    legend(lgd)
    grid on

    %mef(p) = M_EF;
    %g(p) = G;
    figure(2)
    plot(G,M_EF,'g^')
    plot(G,EFavg,'r*')
    ylim([0 1])
    legend('Max Fitness','Average Fitness')
    nme = 'final_population';
    save(nme,'pop_f','G','EF','M_EF');
    G = G+1;
    %p = p+1;
end

saveas(figure(1),'Path','png');
saveas(figure(2),'Fitness','png');
nme = strcat('final_population_G',string(G-1));
save(nme,'pop_f','G','EF','M_EF');

%% Functions

function k = LM(x,y)
    %Determines which directions one can go to stay within a 10x10 grid,
    %given a position on that grid
    if x == 0
        if y == 0
            k = [1 2];
        elseif y == 9
            k = [2 3];
        else
            k = [1 2 3];
        end
    elseif x == 9
        if y == 0
            k = [1 4];
        elseif y == 9
            k = [3 4];
        else
            k = [1 3 4];
        end
    elseif y == 0
        k = [1 2 4];
    elseif y == 9
        k = [2 3 4];
    else
        k = [1 2 3 4];
    end
end