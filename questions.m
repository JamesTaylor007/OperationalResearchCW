%Here, Run this one. This file contains all the different types of questions you might need%
clear
clc
%Our original problem
A = [2, 7, 1, 0, 0, 1;
    5, 8, 0, 2, 0, 0;
    1, 1, 0, 0, 0, 1;
    1, 0, 1, 0, 1, 1];

b = [30, 70, 20, 41];

c = [7, 2, 3, 1, 1, 1];

ineq = [0,-1,1,-1];

minmax = 1;

simplexMethodMatrix(A, b, c, ineq, minmax);


%Another problem, it's Q2, a) in the examples sheet 2
% A2 = [1, 2, -2, 4;
%     2, -1, 1, 2;
%     4, -2, 1, -1];
% 
% b2 = [40; 8; 10];
% c2 = [2, 1, -3, 5];
% ineq2 = [-1, -1, -1];
% minmax2 = 1;
% 
% simplexMethodMatrix(A2, b2, c2, ineq2, minmax2);
    
%Minimising problem
% A3 = [3, 1;
%     4, 3;
%     1, 2];
% 
% b3 = [3; 6; 4];
% c3 = [4, 1];
% ineq3 = [0, 1, -1];
% minmax3 = -1;
% 
% simplexMethodMatrix(A3, b3, c3, ineq3, minmax3);

%Unbounded problem
% A4 = [1, 2, -2, 4;
%     2, -1, 1, 2;
%     4, -2, 2, 4];
% 
% b4 = [40; 8; 10];
% c4 = [2, 1, -3, 5];
% ineq4 = [-1, -1, -1];
% minmax4 = 1;
% 
% simplexMethodMatrix(A4, b4, c4, ineq4, minmax4);

%Infeasible solution
% A5 = [2, 1;
%     3, 4];
% 
% b5 = [2; 12];
% c5 = [3, 2];
% ineq5 = [-1, 1];
% minmax5 = 1;
% 
% simplexMethodMatrix(A5, b5, c5, ineq5, minmax5);
