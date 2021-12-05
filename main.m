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

%Coefficients of z
z = [7+3*M,2+8*M,3+M,1,1,1+2*M,0,0,-M,0,0,-50*M];

simplexMethodMatrix(A,b,z);

%This method would return a vector x containing the solutions x1, x2,
%x3...

%maybe add variables of the BFS to be able to pick the right columns and
%coefficients
function[x] = simplexMethodMatrix(matrix, vector, z_coefficients)

  %xb vector, same as b for now, will change later
  xb = vector;
  cb = sym(zeros(4, 1));
  
  %we should be able to calculate this matrix
  B = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
  B_P_columns = [7,8,10,11];
  
  columns = uint32(1):uint32(length(matrix));
  for i = 1:length(B_P_columns)
    columns = setdiff(columns,B_P_columns(i));
  end
  
  %ENTER STEP 3, OPTIMALITY
  for iteration=1:4
  disp("ITERATION"+iteration)
  %this for loop forms the matrix of those columns in A and coefficients in
  %z
  A_columns = [];
  z_columns = [];
  for i = 1:length(columns)
    A_columns = [A_columns matrix(:,columns(i))];
    z_columns = [z_columns z_coefficients(columns(i))];
  end
  
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
  optimality_vector
  disp("entering P"+entering_col_position);
  
  %ENTER STEP 4, FEASIBILITY
  
  %We use this vector to find which column in B leaves 
  BP_vector = inv(B)*matrix(:,entering_col_position);
  xb_BP_vector = [];
  %Now we find the minimum
  leaving_position_in_B = 1;
  for i = 1:length(BP_vector)
      xb_BP_vector(i) = xb(i)/BP_vector(i);
  end
  xb_BP_vector  
  min_val = max(xb_BP_vector(~isinf(xb_BP_vector)));
  leaving_position_in_B = find(xb_BP_vector==min_val);
  for i=1:length(xb_BP_vector)
      if xb_BP_vector(i)<min_val && 0<= xb_BP_vector(i)
          min_val =  xb_BP_vector(i);
          leaving_position_in_B = i;
      end
  end
leaving_position_in_B
  disp("leaving P"+B_P_columns(leaving_position_in_B));
  
  %ENTER STEP 5, Rewrite values and find z value
  
  %Rewrite B
  B_P_columns(leaving_position_in_B) = entering_col_position
  B_P_columns = sort(B_P_columns)
  temp_position = find(B_P_columns==entering_col_position);
  B = [matrix(:,B_P_columns(1)) matrix(:,B_P_columns(2)) matrix(:,B_P_columns(3)) matrix(:,B_P_columns(4))]
  cb = [z_coefficients(B_P_columns(1)); z_coefficients(B_P_columns(2)); z_coefficients(B_P_columns(3)); z_coefficients(B_P_columns(4))]
  xb = inv(B)*vector;
  %rewrite cb
  cb(temp_position,1) = z_coefficients(entering_col_position)
  current_z = cb.' * xb(1:(length(xb))) + z_coefficients(end)
  end 
end

