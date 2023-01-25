clear
close all
clc

%all time units in minutes
time = 0;
max_time = 8*60;

%total trip time in minutes
av_bike_time = 39;
av_car_time= 44;
av_truck_time = 50;
time_std_dev = 4;

%gallons per minute
gas_bike = 1.7/60;
gas_car = 2.2/60;
gas_truck = 3.2/60;
gas_std_dev = .005;

%capacity per trip
cap_bike = 15;
cap_car = 35;
cap_truck = 50;
cap_std_dev = 4;

%since we don't have actual data i had the code generate it in the
%algorithm as a normal distibution for each variable and they have the same std
%dev to make the process simple but we could change that if we want to
%individualize the variation

%intializing all the "output" variables
raw_data = zeros(1,9);
fitness_values = zeros(1,3);
s_star = [0];
s_w = ["Start"];
min_data = zeros(1,3);


while time < 8*60
    
    %generating values for the fitness function to calculate for
    bike_t = av_bike_time + randn*time_std_dev;
    bike_g = gas_bike + randn*gas_std_dev;
    bike_c = cap_bike + randn*cap_std_dev;

    car_t = av_car_time + randn*time_std_dev;
    car_g = gas_car + randn*gas_std_dev;
    car_c = cap_car + randn*cap_std_dev;

    truck_t = av_truck_time + randn*time_std_dev;
    truck_g = gas_truck + randn*gas_std_dev;
    truck_c = cap_truck + randn*cap_std_dev;
    
    %saving all raw data output points
    current_values = [bike_t,bike_g,bike_c,car_t,car_g,car_c,truck_t,truck_g,truck_c]
    data_hold = raw_data;
    raw_data = [data_hold;current_values]
    
    
    %calculating the gas burned/packages delivered
    bike = (bike_t*bike_g) / bike_c;
    car = (car_t*car_g) / car_c;
    truck = (truck_t*truck_g) / truck_c;
    
    %fitness values matrix saves all the gas/packages values
    current_fit = [bike,car,truck];
    fit_hold = fitness_values;
    fitness_values = [fit_hold;current_fit]
    
    %finding which vehicle minimizes the gas/packages value
    [s,index] = min(current_fit);
    
    if index == 1
        
        %saves the fitness value and index of the lowest value
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + bike_t;
        
        %saves the vehicle type
        s_w_h = s_w;
        s_w = [s_w_h;"Bike"];
       
        %saves the raw data values for time, gas burned, packages for opt
        min_hold = min_data;
        min_data = [min_hold;current_values(1,1),current_values(1,2),current_values(1,3)];
     
        
    elseif index == 2
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + car_t;
        
        s_w_h = s_w;
        s_w = [s_w_h;"Car"];
        
        min_hold = min_data;
        min_data = [min_hold;current_values(1,4),current_values(1,5),current_values(1,6)];
        
        
    else
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + truck_t;
        
        s_w_h = s_w;
        s_w = [s_w_h;"Truck"];
        
        min_hold = min_data;
        min_data = [min_hold;current_values(1,7),current_values(1,8),current_values(1,9)];
    
        
    end
        
end

disp(s_w)

%gas rates
%gas = 3.55; %per gallon
%em = 8887; %gC02/gal of gasoline



