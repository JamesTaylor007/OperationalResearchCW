%Write down the z equation and the constrains

optimising_function = z == 7*x1 + 2*x2 + 3*x3 + x4 + x5 + x6;

%constrains written in canonical form, we don't use them in the code, just
%for clarity
constrain1 = 2*x1 + 7*x2 + x3 + x6 + x7 == 30; 
constrain2 =  5*x1 + 8*x2 + 2*x4 + x8 == 70;
constrain3 = x1 + x2 + x6 - x9 + x10 == 20; 
constrain4 = x1 + x3 + x5 + x6 + x11 == 41;
%This last constrain is to make sure all values are non-negative
constrain5 = x1 * x2 * x3 * x4 * x5 * x6 >= 0;

%We define matrix A with the coefficients of the constrains
A = [2, 7, 1, 0, 0, 1, 1, 0, 0, 0, 0;
     5, 8, 0, 2, 0, 0, 0, 1, 0, 0, 0;
     1, 1, 0, 0, 0, 1, 0, 0, -1, 1, 0;
     1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1];


 
 
 
 
%Method just to test Matlab
disp(method(2));

function[y] = method(x)
    y = 2*x;
end

