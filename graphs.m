%Write down the z equation and the constrains
syms z x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 M;
optimising_function = z == 7*x1 + 2*x2 + 3*x3 + x4 + x5 + x6;

%constrains written in canonical form, we don't use them in the code, just
%for clarity
constrain1 = sym('constrain1_%d',[1 100]);
for i=1:length(constrain1)
    constrain1(i) = 2*x1 + 7*x2 + x3 + x6 + x7 == i;
end

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
b_1 = [1:1:100];
b_2 = 70;
b_3 = 20;
b_4 = 41;


%Coefficients of z
z = [7+3*M,2+8*M,3+M,1,1,1+2*M,0,0,-M,0,0,-50*M];

answer = zeros(100);
disp([b_1(5);b_2;b_3;b_4])
for i=1:length(answer)
    simplexMethodMatrix(A,[i;b_2;b_3;b_4],z)
end


%This method would return a vector x containing the solutions x1, x2,
%x3...

%maybe add variables of the BFS to be able to pick the right columns and
%coefficients
function z_solution = simplexMethodMatrix(matrix, vector, z_coefficients)

  z_solution = 0;
  x_values = [];
  x_solutions = [];

  %xb vector, same as b for now, will change later
  xb = vector;
  cb = sym(zeros(4, 1));
  
  %B matrix
  B = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
  %This variable will be used to know which P's are currently in B
  B_P_columns = [7,8,10,11];

  %The columns variable will be used to know which columns are in A and
  %which columns are in B
  columns = uint32(1):uint32(length(matrix));
  for i = 1:length(B_P_columns)
    columns = setdiff(columns,B_P_columns(i));
  end
  
  %This variable will be used to check when the optimum has been found,
  %when it is positive, the while loop will stop
  all_positive = -1;
  while all_positive<0
 
  %this for loop forms the matrix of those columns in A
  %and coefficients in z
  A_columns = [];
  z_columns = [];
  for i = 1:length(columns)
    A_columns = [A_columns matrix(:,columns(i))];
    z_columns = [z_columns z_coefficients(columns(i))];
  end
  
  %OPTIMALITY STEP
  
  optimality_vector = cb.' * inv(B) * A_columns - z_columns;

  %we find out the most negative value and save the position so we know
  %which vector is going to enter B
  entering_col_position = 0;
  min_value = 0;
  format long;
  for i = 1:length(optimality_vector)
      temp = optimality_vector(i);
      syms M;
      current_value = double(subs(temp,M,1000000));
      if current_value<min_value
          entering_col_position = columns(i);
          min_value = current_value;
      end
  end
  
  %If there are no negative values, a solution has been found
  all_positive = min_value;
  if all_positive<0
    
    %FEASIBILITY STEP
    
    %We use this vector to find which column in B leaves
    BP_vector = inv(B)*matrix(:,entering_col_position);
    xb_BP_vector = [];
    %Now we find the minimum
    leaving_position_in_B = 1;
    for i = 1:length(BP_vector)
        xb_BP_vector(i) = xb(i)/BP_vector(i);
    end
    min_val = max(xb_BP_vector(~isinf(xb_BP_vector)));
    leaving_position_in_B = find(xb_BP_vector==min_val);
    for i=1:length(xb_BP_vector)
        if xb_BP_vector(i)<min_val && 0<= xb_BP_vector(i)
              min_val =  xb_BP_vector(i);
              leaving_position_in_B = i;
        end
    end

      %REWRITE VALUES AND FIND Z
      
      %Rewrite B
      B_P_columns(leaving_position_in_B) = entering_col_position;
      B_P_columns = sort(B_P_columns);
      B = [matrix(:,B_P_columns(1)) matrix(:,B_P_columns(2)) matrix(:,B_P_columns(3)) matrix(:,B_P_columns(4))];
      cb = [z_coefficients(B_P_columns(1)); z_coefficients(B_P_columns(2)); z_coefficients(B_P_columns(3)); z_coefficients(B_P_columns(4))];
      xb = inv(B)*vector;
      %rewrite cb
      current_z = cb.' * xb(1:(length(xb))) + z_coefficients(end);
  end
  end
  x_solutions = xb;
  z_solution = current_z;
  x_values = B_P_columns;
  
end





