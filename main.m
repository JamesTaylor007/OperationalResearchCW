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




function[valid] = checkIfValidInput(constrainsMatrix, b_values, z_coefficients, inequalities)
  %The valid variable stores whether the input is valid or not, if it's 1 
  %it's valid and if it's -1 it's not
  valid = 1;
  
  %First we check for redundant equations
  n_constrains = length(constrainsMatrix(:,1));
  n_var_in_A = length(constrainsMatrix(1,:));
  n_b_values = length(b_values);
  n_ineq = length(inequalities);
  n_coefficients = length(z_coefficients);
  
  if n_ineq~=n_b_values || n_ineq~=n_constrains
     valid = -1; 
  end
  
  if n_coefficients~=n_var_in_A
     valid = -1; 
  end
  
  if valid==-1
     disp("Size of inputs doesn't match,")
     disp("please check again before running the program");
     return;
  end
  
  for i=1:n_b_values
     for j=1:n_b_values
        if i~=j
          scalar = 0;
          temp_vector_1 = constrainsMatrix(i,:);
          temp_vector_2 = constrainsMatrix(j,:);
          for k=1:length(temp_vector_1)
              if temp_vector_1(k)~=0
                  scalar = temp_vector_1(k)/temp_vector_2(k);
              end
          end
          valid = -1;
          for k=1:length(temp_vector_1)
              if scalar*temp_vector_2(k)~=temp_vector_1(k)
                  valid = 1;
              end
          end 
          if valid == -1 
              disp("There's two redundant equations,");
              disp("Please check your input");
              return
          
          end
        end
     end
  end
end

function[constrainsMatrix, B, b_values, b_columns, z_coefficients] = turnToCanonicalForm(constrainsMatrix, b_values, z_coefficients, inequalities, minmax)
  initial_length = length(constrainsMatrix(1,:));
  artificial_var = [];
  n_of_vars = initial_length;

  %If b_vector is not in column form it turns it into one
   if length(b_values(1,:))~=1
      b_values = b_values.';
   end
 
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
    if minmax==1
        temp_row = constrainsMatrix(temp_number,:) * M;
    else
        constrainsMatrix(temp_number,:);
        temp_row = -constrainsMatrix(temp_number,:) * M;
    end
    for j=1:length(b_columns)
        temp_row(b_columns(j)) = 0;
    end

    z_coefficients = z_coefficients + temp_row;
    constant = constant - b_values(temp_number);
  end
  if minmax==1
      z_coefficients = [z_coefficients constant*M];
  else
      z_coefficients = [z_coefficients -constant*M];
  end
  
  disp("The problem in canonical form becomes:");
  disp(constrainsMatrix);
  
  disp("With the z coefficients being:");
  disp(z_coefficients);
  
end

function[] = simplexMethodMatrix(constrainsMatrix, b_values, z_coefficients, inequalities, minmax)

  z_solution = 0;
  x_values = [];
  x_solutions = [];
  current_z = 0;
  past_zs = [];
  
  valid = checkIfValidInput(constrainsMatrix, b_values, z_coefficients, inequalities);
  if valid==-1
     return;
  end
  [matrix, B, b_values, b_columns, z_coefficients] = turnToCanonicalForm(constrainsMatrix, b_values, z_coefficients, inequalities, minmax);
  
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
  count = -1;
  while all_positive<0
      count = count + 1;
      disp("current iteration: "+count);
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
  min_max_value = 0;
  format long;
  big_M = 1000;
  for i = 1:length(optimality_vector)
      temp = optimality_vector(i);
      syms M;
      current_value = double(subs(temp,M,big_M));
      if minmax==1
          if current_value<min_max_value
              entering_col_position = columns(i);
              min_max_value = current_value;
          end
      else
          if current_value>min_max_value
              entering_col_position = columns(i);
              min_max_value = current_value;
          end
      end
  end

  %If there are no negative values, a solution has been found
  all_positive = min_max_value;
  if minmax==-1
      all_positive = -all_positive;
  end

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
    if minmax==-1
        min_val = max(xb_BP_vector(~isinf(xb_BP_vector)));
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
      
      %If the last z is equal to the current z, it means that we've found 
      %the optimum
      if any(past_zs(:) == current_z(end))
          all_positive = 1;
      else
          if any(past_zs(:) == current_z)
            disp('Problem is cyclic');
            all_positive = 1;
          else
            past_zs = [past_zs current_z];
          end
      end

  end
  end
  
  x_solutions = xb;
  z_solution = current_z;

  %Infeasible solutions
  if has(z_solution, M)
     disp(' '); 
     disp('This problem is unfeasible'); 
  else
      
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
end

