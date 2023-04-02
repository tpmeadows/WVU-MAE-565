function u_star = test_defuzz(irm_fuzz)

%trapezoid base sizes
O_LB = 24;
O_SB = 8;

I_LB = 24;
I_SB = 8;

LB_mat = [I_LB,I_LB,I_LB,O_LB,O_LB;I_LB,I_LB,I_LB,I_LB,O_LB;I_LB,I_LB,I_LB,I_LB,I_LB;O_LB,I_LB,I_LB,I_LB,O_LB;O_LB,O_LB,I_LB,I_LB,I_LB];
SB_mat = [I_SB,I_SB,I_SB,O_SB,O_SB;I_SB,I_SB,I_SB,I_SB,O_SB;I_SB,I_SB,I_SB,I_SB,I_SB;O_SB,I_SB,I_SB,I_SB,O_SB;O_SB,O_SB,I_SB,I_SB,I_SB];

%area calulcations
area_mat = zeros(5,5);

a = 1;
b = 1;

while a < 6
    while b < 6
        area_mat(a,b) = LB_mat(a,b)*irm_fuzz(a,b) - (((irm_fuzz(a,b))^2)/2)*(LB_mat(a,b) - SB_mat(a,b));
        b = b + 1;
        
    end
    b = 1;
    a = a + 1;
end
disp(area_mat)

%centroids
CLN = -24;
CN = -12;
CZ = 0;
CP = 12;
CLP = 24;

cent_mat = [CZ,CP,CP,CLP,CLP;CN,CZ,CP,CP,CLP;CN,CN,CZ,CP,CP;CLN,CN,CN,CZ,CP;CLN,CLN,CN,CN,CZ];

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
disp(numerator)

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
disp(denom)

u_star = numerator/denom

end