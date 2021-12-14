%clear your work space and console so that there is no clash
clear
clc
%Write down the z equation and the constrains
syms z x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 M;

A = [2, 7, 1, 0, 0, 1;
    5, 8, 0, 2, 0, 0;
    1, 1, 0, 0, 0, 1;
    1, 0, 1, 0, 1, 1];

b = [30; 70; 20; 41];

coefficients = [7, 2, 3, 1, 1, 1];

inequalities = [0,-1,1,-1];

minmax = 1;



%----------------%
% --Contraint 1--%
%----------------%

bValues = [20:1:33];

zSolutions_C1 = zeros(1,length(bValues));
for i=1:length(zSolutions_C1)
    zSolutions_C1(i) = simplexMethodMatrix(A, [bValues(i); 70; 20; 41], coefficients, inequalities, minmax);
end
disp(zSolutions_C1)
figure(1)
plot(bValues,zSolutions_C1)
title("Z values changing in response to contraint 1")
xlabel('Constraint value') 
ylabel('z value') 


%----------------%
% --Contraint 2--%
%----------------%

bValues1 = [52:1:100];

zSolutions_C2 = zeros(1,length(bValues1));
for i=1:length(zSolutions_C2)
    zSolutions_C2(i) = simplexMethodMatrix(A, [30; bValues1(i); 20; 41], coefficients, inequalities, 1);
end
disp(zSolutions_C2)
figure(2)
plot(bValues1,zSolutions_C2)
title("Z values changing in response to contraint 2")
xlabel('Constraint value') 
ylabel('z value') 


%----------------%
% --Contraint 3--%
%----------------%

bValues2 = [18:1:30];

zSolutions_C3 = zeros(1,length(bValues2));
for i=1:length(zSolutions_C3)
    zSolutions_C3(i) = simplexMethodMatrix(A, [30; 70; bValues2(i); 41], coefficients, inequalities, 1);
end
disp(zSolutions_C3)
figure(3)
plot(bValues2,zSolutions_C3)
title("Z values changing in response to contraint 3")
xlabel('Constraint value') 
ylabel('z value') 


%----------------%
% --Contraint 4--%
%----------------%

bValues3 = [21:1:60];

zSolutions_C4 = zeros(1,length(bValues3));
for i=1:length(zSolutions_C4)
    zSolutions_C4(i) = simplexMethodMatrix(A, [30; 70; 20; bValues3(i)], coefficients, inequalities, 1);
end
disp(zSolutions_C4)
figure(4)
plot(bValues3,zSolutions_C4)
title("Z values changing in response to contraint 4")
xlabel('Constraint value') 
ylabel('z value') 







