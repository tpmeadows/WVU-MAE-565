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