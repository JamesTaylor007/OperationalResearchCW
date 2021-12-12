% We define matrix A with the coefficients of the constrains
% A = [2, 7, 1, 0, 0, 1, 1, 0, 0, 0, 0;
%      5, 8, 0, 2, 0, 0, 0, 1, 0, 0, 0;
%      1, 1, 0, 0, 0, 1, 0, 0, -1, 1, 0;
%      1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1];
% 
% b vector 
% b = [30; 70; 20; 41];
% 
% the columns of A which we take to form B
% cols = [7,8,10,11];
% 
% Coefficients of z
% z = [7+3*M,2+8*M,3+M,1,1,1+2*M,0,0,-M,0,0,-50*M];

A = [2, 7, 1, 0, 0, 1;
    5, 8, 0, 2, 0, 0;
    1, 1, 0, 0, 0, 1;
    1, 0, 1, 0, 1, 1];

b = [30; 70; 20; 41];

coefficients = [7, 2, 3, 1, 1, 1];

inequalities = [0,-1,1,-1];

minmax = 1;

%Another problem

A2 = [2, 1, 1;
    1, 1, 0;
    1, 2, 0];

b2 = [20; 18; 12];
coefficients2 = [5, 4, 0];
inequalities2 = [-1,-1,1];

simplexMethodMatrix(A, b, coefficients, inequalities, 1);
%simplexMethodMatrix(A2,b2,coefficients2,inequalities2,1);

function[constrainsMatrix, B, b_columns, z_coefficients] = turnToCanonicalForm(constrainsMatrix, b_values, z_coefficients, inequalities)
  initial_length = length(constrainsMatrix);
  artificial_var = [];
  n_of_vars = initial_length;
  for i=1:length(inequalities)
      n_of_vars = n_of_vars + 1;
      %If less or equal sign
      if inequalities(i)==-1
          %Add slack variable
          col_to_add = zeros(length(constrainsMatrix(:,1)), 1);
          col_to_add(i) = 1;
          constrainsMatrix = [constrainsMatrix col_to_add];
          z_coefficients = [z_coefficients 0];
      end
      if inequalities(i)==1
          n_of_vars = n_of_vars + 1;
          %Add slack variable and artificial
          col_to_add = zeros(length(constrainsMatrix(:,1)), 1);
          col_to_add2 = zeros(length(constrainsMatrix(:,1)), 1);
          col_to_add(i) = -1;
          col_to_add2(i) = 1;
          constrainsMatrix = [constrainsMatrix col_to_add];
          constrainsMatrix = [constrainsMatrix col_to_add2];
          artificial_var = [artificial_var n_of_vars];
          z_coefficients = [z_coefficients 0 0];
      end
      if inequalities(i)==0
          %Add artificial variable
          col_to_add = zeros(length(constrainsMatrix(:,1)), 1);
          col_to_add(i) = 1;
          constrainsMatrix = [constrainsMatrix col_to_add];
          artificial_var = [artificial_var n_of_vars];
          z_coefficients = [z_coefficients 0];
      end
  end

  B = [];
  b_columns = [];
  %Now we find the BFS and columns in A which we will use in B
  count = 1;
  for i=(initial_length+1):length(constrainsMatrix)
      col_to_add = zeros(length(constrainsMatrix(:,1)), 1);
      col_to_add(count) = 1;
      if constrainsMatrix(:,i)==col_to_add
        B = [B col_to_add];
        b_columns = [b_columns i];
        count = count + 1;
      end
  end

  constant = 0;
  syms M;
  %Transform z equation
  for i=1:length(artificial_var)
    temp_col = constrainsMatrix(:,artificial_var(i));
    temp_number = find(temp_col==1);
    temp_row = sym(zeros(length(constrainsMatrix), 1));
    temp_row = constrainsMatrix(temp_number,:) * M;
    for j=1:length(b_columns)
        temp_row(b_columns(j)) = 0;
    end
    z_coefficients = z_coefficients + temp_row;
    constant = constant - b_values(temp_number);
  end

  z_coefficients = [z_coefficients constant*M];
  
end

%This method would return a vector x containing the solutions x1, x2,
%x3...

%maybe add variables of the BFS to be able to pick the right columns and
%coefficients
function[] = simplexMethodMatrix(constrainsMatrix, b_values, z_coefficients, inequalities, minmax)

  z_solution = 0;
  x_values = [];
  x_solutions = [];
  
  [matrix, B, b_columns, z_coefficients] = turnToCanonicalForm(constrainsMatrix, b_values, z_coefficients, inequalities);
  
  %This variable will be used to know which P's are currently in B
  B_P_columns = b_columns;

  %xb vector, same as b for now, will change later
  xb = b_values;
  cb = sym(zeros(length(xb), 1));
  for i=1:length(B_P_columns)
    cb(i,1) = z_coefficients(B_P_columns(i));
  end
  

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
  big_M = 0;
  if minmax==1
     big_M = 1000000;
  else
     big_M = -1000000;
  end
  for i = 1:length(optimality_vector)
      temp = optimality_vector(i);
      syms M;
      current_value = double(subs(temp,M,big_M));
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
      for i=1:length(B_P_columns)
        B(:,i) = matrix(:,B_P_columns(i));
        cb(i,1) = z_coefficients(B_P_columns(i));
      end
      
      xb = inv(B)*b_values;
      
      current_z = cb.' * xb(1:(length(xb))) + z_coefficients(end);
  end
  end
  x_solutions = xb;
  z_solution = current_z;
  x_values = B_P_columns;
  disp(" ")
  disp("solution found:")
  disp(" ")
  disp("z value:")
  disp(z_solution);
  disp("x values:")
  for i=1:length(x_solutions)
      disp("x"+x_values(i)+" = "+xb(i))
  end
end

