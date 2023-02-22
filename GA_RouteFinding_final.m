%Route-finding GA
%Last edited: 2-21-23 @ 8:07pm

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
d = 1; %distance between points, miles

n = 100; %pop size
moves = [0 1;1 0;0 -1;-1 0];
m = round(0.10*n); %# of mutating individuals
c = round(0.75*n); %# of reproducing individuals
if rem(c,2) ~= 0
    c = c + 1;
end

t = 12; %work day, hours
L0 = t*v+2; %max # of stops
chk0 = 100; %# of visited houses - origin
pop(1:n) = struct('path',zeros(L0,2),'turns',zeros(L0,1),'distance',0,'delivered',0,'stock',zeros(L0,1));

%Criteria
deliver0 = 120; %# of deliveries
d0 = 541;
criteria = [d0,deliver0];

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
    tr = t; %max time, hours
    chk = chk0;
    delivered = 0;
    stops = [0 0];
    carry = cap;
    a = []; %storage variable
    direction = zeros(L0,1);
    stock = zeros(L0,1);
    total_distance = 0;
    j = 1;
    while tr > 0
       %for directions, N=1, E=2, S=3, W=4
       %k vector represents the set of legal moves
       x = stops(j,1); y = stops(j,2);
       if carry == 0
            while (x ~= 0 || y ~= 0) && tr > 0
                if x ~= 0
                    direction(j) = 4;
                    x = x - 1;
                    j = j+1;
                    stops(j,:) = [x y];
                elseif y ~= 0
                    direction(j) = 3;
                    y = y - 1;
                    j = j+1;
                    stops(j,:) = [x y];
                end
                 %d = sqrt((stops(j-1,1)-x)^2 + (stops(j-1,2)-y)^2);
                 total_distance = total_distance + d;
                 tr = tr - d/v;
            end
            carry = cap;
       else
           if chk == 0
               chk = chk0;
               a = [a;stops];
               stops = stops(end,:);
           end
           
           if x == 0 && y == 0
               carry = cap;
           end

           k = LM(x,y);
           l = length(k);
           direction(j) = k(randi(l));
           new_stop = [x,y] + moves(direction(j),:);
           if isempty(find(stops(:,1) == new_stop(1) & stops(:,2) == new_stop(2),1)) == true ...
                   && (x ~= 0 || y ~= 0)
               chk = chk - 1;
               carry = carry - 1;
               delivered = delivered + 1;
           end
           stops(j+1,:) = new_stop;
           stock(j) = carry;
           j = j+1;
           %d = sqrt((stops(j,1)-x)^2 + (stops(j,2)-y)^2);
           total_distance = total_distance + d;
           tr = tr - d/v;
       end
    end
    pop(i).path = stops;%[a;stops];
    k = LM(stops(j,1),stops(j,2));
    direction(j) = k(randi(length(k)));
    pop(i).turns = direction;
    pop(i).distance = total_distance;
    pop(i).delivered = delivered;
    pop(i).stock = stock;
    EF(i,1) = Fitness([total_distance,delivered],criteria);
end

M_EF = max(EF);
ind_max = find(EF == M_EF);
ind_max = ind_max(randi(length(ind_max)));
G = 1;
p = 1;
while M_EF <= 0.9

%%%%%%%%% Selection (Roulette Wheel, Elitist) %%%%%%%%%
    
    TF = sum(EF);
    q = zeros(n,1);
    for i = 1:n
        q(i) = sum(EF(1:i))/TF;
    end
    pop0(1:n) = struct('path',zeros(L0,2),'turns',zeros(L0,1),'distance',0,'delivered',0,'stock',zeros(L0,1));
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
        if ind_m(i) ~= ind_max
            excl = find(~pop_m(i).stock);
            pop_m(i).path(excl,:) = [];
            pop_m(i).turns(excl) = [];
            pop_m(i).stock(excl) = [];

            Lm = length(pop_m(i).turns);
            r = randi(Lm);
            pos_m = pop_m(i).path(r,:);
            turn_m = pop_m(i).turns(r);

            k = LM(pos_m(1),pos_m(2));
            ind = find(k ~= turn_m);
            l = length(ind);
            mv_m = k(ind(randi(l)));
            pop_m(i).turns(r) = mv_m;
            j = r;
            while j ~= Lm && (pop_m(i).path(j+1,1) ~= 0 || pop_m(i).path(j+1,2) ~= 0)
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
                j = j+1;
            end
            
%             if s ~= L0
%                 r = incl(find(incl > s,1));
%                 path_str = pop_m(i).path(r:end,:);
%                 pop_m(i).path(s+1:end,:) = [];
%                 turn_str = pop_m(i).turns(r:end);
%                 pop_m(i).turns(s+1:end) = [];
% 
%                 direction = [];
%                 stops = [0 0];
%                 x = pos_m(1); y = pos_m(2);
%                 p = 1;
%                 while x ~= 0 || y ~= 0
%                     if x ~= 0
%                         direction(p,1) = 4;
%                         x = x - 1;
%                         stops(p,:) = [x y];
%                         p = p+1;
%                     elseif y ~= 0
%                         direction(p,1) = 3;
%                         y = y - 1;                        
%                         stops(p,:) = [x y];
%                         p = p+1;
%                     end
%                 end
%                 pop_m(i).path = [pop_m(i).path;stops;path_str];
%                 pop_m(i).turns = [pop_m(i).turns;direction;turn_str];
%                 Lmf = length(pop_m(i).turns);
%                 if Lmf > L0
%                     pop_m(i).path = pop_m(i).path(1:L0,:);
%                     pop_m(i).turns = pop_m(i).turns(1:L0,:);
%                 elseif Lmf < L0
%                     pos_m = pop_m(i).path(end,:);
%                     turn_m = pop_m(i).turns(end);
% 
%                     k = LM(pos_m(1),pos_m(2));
%                     ind = find(k ~= turn_m);
%                     l = length(ind);
%                     mv_m = k(ind(randi(l)));
%                     pop_m(i).turns(end) = mv_m;
%                     for j = Lmf:L0-1
%                         b = pos_m + moves(mv_m,:);
%                         if b(1) > 9 || b(1) < 0 || b(2) > 9 || b(2) < 0
%                             k = LM(pos_m(1),pos_m(2));
%                             ind = find(k ~= turn_m);
%                             l = length(ind);
%                             mv_m = k(ind(randi(l)));
%                             pop_m(i).turns(j) = mv_m;
%                             b = pos_m + moves(mv_m,:);
%                         end
%                         pop_m(i).path(j+1,:) = b;
%                         pos_m = b;
%                         k = LM(pos_m(1),pos_m(2));
%                         ind = find(k ~= turn_m);
%                         l = length(ind);
%                         mv_m = k(ind(randi(l)));
%                         pop_m(i).turns(j+1) = mv_m;
%                     end
%                 end
%             end
            pop0(ind_m(i)) = pop_m(i);
        end
    end
    %end
    
    %%%%%%%%%    Crossover    %%%%%%%%%
    ind_c = randi(n,c,1); %indices of reproducing individuals
    if isempty(find(ind_c == ind_max,1)) == false
        l_max = find(ind_c == ind_max);
        ind_c(l_max) = [];
        for i = 1:length(l_max)
            lic = length(ind_c);
            ind_c(randi(lic)) = [];
        end
    end
    lic = length(ind_c);
    pop_c = pop0(ind_c);

    for i = 1:lic/2
        i1 = ind_c(randi(lic/2));
        i2 = ind_c(randi([round(lic/2)+1 lic]));
        %if i1 ~= ind_max && i2 ~= ind_max
            pos1 = pop0(i1).path;
            turn1 = pop0(i1).turns;
            stock1 = pop0(i1).stock;
            ign1 = find(~stock1);
            pos1(ign1,:) = []; turn1(ign1) = []; stock1(ign1) = [];
            l1 = length(turn1);

            pos2 = pop0(i2).path;
            turn2 = pop0(i2).turns;
            stock2 = pop0(i2).stock;
            ign2 = find(~stock2);
            pos2(ign2,:) = []; turn2(ign2) = []; stock2(ign2) = [];
            l2 = length(turn2);

            if l1 > l2
                i_crs = find(~(pos1(1:l2,1)-pos2(:,1)) & ~(pos1(1:l2,2)-pos2(:,2)));
                intvl = [round(l2*0.1) round(l2*0.9)]; %bounds for crossover point
            elseif l1 <= l2
                i_crs = find(~(pos1(:,1)-pos2(1:l1,1)) & ~(pos1(:,2)-pos2(1:l1,2)));
                intvl = [round(l1*0.1) round(l1*0.9)]; %bounds for crossover point
            end
            i_crs = i_crs(i_crs > intvl(1) & i_crs < intvl(2));
            lc = length(i_crs);
            if lc ~= 0
                crosspt = i_crs(randi(lc));
                a = pos1; at = turn1; as = stock1;
                b = pos2; bt = turn2; bs = stock2;
                cc = a(crosspt:end,:); ct = at(crosspt:end,:);
                cs = as(crosspt:end);
                a(crosspt:end,:) = []; a = [a;b(crosspt:end,:)];
                at(crosspt:end,:) = []; at = [at;bt(crosspt:end,:)];
                as(crosspt:end,:) = []; as = [as;bs(crosspt:end,:)];

                b(crosspt:end,:) = []; b = [b;cc];
                bt(crosspt:end,:) = []; bt = [bt;ct];
                bs(crosspt:end,:) = []; bs = [bs;cs];

                pos1 = a; turn1 = at; stock1 = as;
                pos2 = b; turn2 = bt; stock2 = bs;
            end
            pop_c(i).path = pos1;
            pop_c(i).turns = turn1;
            pop_c(i).stock = stock1;
            pop_c(i+c/2-1).path = pos2;
            pop_c(i+c/2-1).turns = turn2;
            pop_c(i+c/2-1).stock = stock2;
%         else
%             pop_c(i) = pop0(i1);
%             pop_c(i+c/2-1) = pop0(i2);
%         end
    end
    pop0(ind_c) = pop_c(1:length(ind_c));

%%%%%%%%%    Fitness Evaluation     %%%%%%%%%
    for i = 1:n
        chk = chk0;
        carry = cap;
        delivered = 0;
        total_distance = 0;
        stops = [];
        x = []; y = [];
        c_stops = 1;
        stock = zeros(L0,1);
        a = 1;
        b = 0;
        tr = t + d/v;
        j = 1;
        while tr > 0
           
            homes = find(pop0(i).path(:,1) == 0 & pop0(i).path(:,2) == 0);
            if carry == 0
                x = stops(c_stops,1); y = stops(c_stops,2);
                cutoff = homes(find(homes > j-1,1));
                if isempty(cutoff) == false
                    path_str = pop0(i).path(cutoff:end,:);
                    turn_str = pop0(i).turns(cutoff:end);
                end
                pop0(i).path(j:end,:) = [];
                pop0(i).turns(j:end) = [];
                while (x ~= 0 || y ~= 0) && tr > 0
                    if x ~= 0
                        x = x - 1;
                        pop0(i).path(j,:) = [x y];
                        pop0(i).turns(j) = 4;
                        j = j+1;
                    elseif y ~= 0
                        y = y - 1;
                        pop0(i).path(j,:) = [x y];
                        pop0(i).turns(j) = 3;      
                        j = j+1;
                    end
                    %d = sqrt((stops(j-1,1)-x)^2 + (stops(j-1,2)-y)^2);
                    total_distance = total_distance + d;
                    tr = tr - d/v;
                end
                pop0(i).path = [pop0(i).path;path_str(2:end,:)];
                pop0(i).turns = [pop0(i).turns;turn_str(2:end)];
                carry = cap;
            else
                if j > length(pop0(i).turns)
                    k = LM(x,y);
                    l = length(k);
                    direction = k(randi(l));
                    pop0(i).path(j,:) = [x,y] + moves(direction,:);
                    pop0(i).turns(j) = direction;
                end
                stops = pop0(i).path(a:j,:);
                c_stops = j-b;
                x = stops(c_stops,1); y = stops(c_stops,2);

                if j ~= 1 %< length(pop0(i).turns)
                    if abs(x - pop0(i).path(j-1,1)) > 1 || abs(y - pop0(i).path(j-1,2)) > 1 ...
                            || (abs(x - pop0(i).path(j-1,1)) == 1 && abs(y - pop0(i).path(j-1,2)) == 1)
                        path_str = pop0(i).path(j:end,:);
                        pop0(i).path(j:end,:) = [];
                        turn_str = pop0(i).turns(j:end);
                        pop0(i).turns(j:end) = [];
                        
                        x0 = pop0(i).path(j-1,1); y0 = pop0(i).path(j-1,2);
                        k = LM(x0,y0);
                        l = length(k);
                        direction = k(randi(l));
                        new_stop = [x0 y0] + moves(direction,:);
                        stops(j,:) = new_stop;
                        pop0(i).turns = [pop0(i).turns;direction;turn_str];
                        pop0(i).path = [pop0(i).path;new_stop;path_str];
                    end
                end       

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

                %d = 1/v*sqrt((stops(j-1,1)-x)^2 + (stops(j-1,2)-y)^2);
                total_distance = total_distance + d;
                tr = tr - d/v;                
                j = j+1;
            end
            stock(j) = carry;
            path_str = []; turn_str = [];
        end

        pop0(i).distance = total_distance;
        pop0(i).delivered = delivered;
        pop0(i).stock = stock;

        if length(pop0(i).turns) > L0 || length(pop0(i).stock) > L0
            pop0(i).path(L0+1:end,:) = [];
            pop0(i).turns(L0+1:end) = [];
            pop0(i).stock(L0+1:end) = [];
        end

        EF(i,1) = Fitness([total_distance,delivered],criteria);
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

function EvFit = Fitness(param,crit)
%Parameter #1 is minimized (EF(d = 200) = 1, EF(d = dist0) = 0)
%Parameter #2 is maximized (EF(deliv = deliv0) = 1, EF(deliv = 0) = 0)
    dist = param(1); 
    deliv = param(2);
    dist0 = crit(1); 
    deliv0 = crit(2);

    EF_dist = (1/(200-dist0))*(dist-dist0);
    EF_deliv = 1 + (1/deliv0)*(deliv-deliv0);
    
    EvFit = 0.3*EF_dist + 0.7*EF_deliv;
end