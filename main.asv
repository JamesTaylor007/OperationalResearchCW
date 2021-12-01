%Write down the z equation and the constrains
syms z x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 M;
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

%b vector 
b = [30; 70; 20; 41];

%Might be useful to define z as a vector of its coefficients
z = [7+3*M,2+8*M,3+M,1,1,1+2*M,0,0,M,0,0,-50*M];

%This method would return a vector x containing the solutions x1, x2,
%x3...

%maybe add variables of the BFS to be able to pick the right columns and
%coefficients
function[x] = simplexMethodMatrix(matrix, vector)

  %xb vector, same as b for now, will change later
  xb = b;
  cb = [0; 0; 0; 0];
  
  %we should be able to calculate this matrix
  B = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]
  
end

