clear
close all
clc

%all time units in minutes
time = 0;
max_time = 8*60;

%total trip time in minutes
av_bike_time = 34;
av_car_time= 39;
av_truck_time = 46;
time_std_dev = 5;

%gallons per minute
gas_bike = 1.7/60;
gas_car = 2.2/60;
gas_truck = 3.4/60;
gas_std_dev = .005;

%capacity per trip
cap_bike = 15;
cap_car = 35;
cap_truck = 50;
cap_std_dev = 4;


raw_data = zeros(1,9);
fitness_values = zeros(1,3);
s_star = [0];
s_w = ["Start"];

while time < 8*60
    
    bike_t = av_bike_time + randn*time_std_dev;
    bike_g = gas_bike + randn*gas_std_dev;
    bike_c = cap_bike + randn*cap_std_dev;

    car_t = av_car_time + randn*time_std_dev;
    car_g = gas_car + randn*gas_std_dev;
    car_c = cap_car + randn*cap_std_dev;

    truck_t = av_truck_time + randn*time_std_dev;
    truck_g = gas_truck + randn*gas_std_dev;
    truck_c = cap_truck + randn*cap_std_dev;

    
    bike = (bike_t*bike_g) / bike_c;
    car = (car_t*car_g) / car_c;
    truck = (truck_t*truck_g) / truck_c;
    
    current_values = [bike_t,bike_g,bike_c,car_t,car_g,car_c,truck_t,truck_g,truck_c]
    data_hold = raw_data;
    raw_data = [data_hold;current_values]
    
    current_fit = [bike,car,truck];
    fit_hold = fitness_values;
    fitness_values = [fit_hold;current_fit]
    
    [s,index] = min(current_fit);
    
    if index == 1
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + bike_t;
        
        s_w_h = s_w;
        s_w = [s_w_h;"Bike"];
       
        
    elseif index == 2
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + car_t;
        
        s_w_h = s_w;
        s_w = [s_w_h;"Car"];
        
    else
        s_hold = s_star;
        s_star = [s_hold;index];
        time = time + truck_t;
        
        s_w_h = s_w;
        s_w = [s_w_h;"Truck"];
        
    end
        
end

disp(s_w)

%gas rates
gas = 3.55; %per gallon
em = 8887; %gC02/gal of gasoline



