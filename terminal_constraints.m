function [F, c, P] = terminal_constraints(system, Q, R, constraints)

% TERMINAL CONSTRAINTS
%
% [F,c,P] = terminal_constraints(system, Q, R, constraints) calculates the
% optimal set of constraints on the terminal state vector of a
% discrete-time state space system. Using the Q and R matrices, this
% function resturns a set of terminal constraints and a matrix for use in
% computing terminal weights
%
% This will allow the lqr controller to go towards the origin without violating 
% the constraints
% 
% Note:
%   * Q matrix is a weighted sum of squares on the states so the size of the Q matrix
%     needs to be the same as the A matrix
%   * R matrix size needs to be that of the amount of inputs into the system
%
%
% Written By Y3843416
%

%% Extracting system matrices

A = system.A;
B = system.B;
C = system.C;
D = system.D;

%% Check Q and R sizes

% Q is a weighted sum of squares on the states so the size of the Q matrix
% needs to be the same as the A matrix. R is a weighted sum of squares of
% the inputs, therefor, the size of R needs to be that of the amount of
% inputs into the system. If these sizes are not correct then an error
% needs to be thrown

if size(Q,1) ~= size(A,1) || size(Q,2) ~= size(A,2)
    msg = ['Error Occured. Q matrix is not the correct size. Must be same dimensions as the number of states: ' num2str(size(A,1))];
    error(msg)
end
if size(R,2) ~= size(B,2) || size(R,1) ~= size(B,2)
    msg = ['Error Occured. R matrix is not the correct size. Must be the same size as the number of inputs: ' num2str(size(B,2))];
    error(msg)
end

%% Discrete LQR

% Computing the K and Symmetric P matrix using discrete LQR
[K,P] = dlqr(A,B,Q,R);

%% Inital single sided constraints

% Creating the initail single sided constraint matrix and vector.
% Using if statements to look through the structure and check if the
% correct field is present. If the field does not exist then the F matrix
% and c vector do not get populated. 

% Initialising c vector and F matrix
c = [];
F = [];

if isfield(constraints, 'umin')
    F = [F;-K];
    c = abs([c;constraints.umin]);
end
if isfield(constraints, 'umax')
    F = [F;K];
    c = abs([c;constraints.umax]);
end
if isfield(constraints, 'xmin') 
    F = [F;-eye(size(A))];
    c = abs([c;constraints.xmin]);
end
if isfield(constraints, 'xmax') 
    F = [F;eye(size(A))];
    c = abs([c;constraints.xmax]);
end
if isfield(constraints, 'ymin')
    F = [F;-C];
    c = abs([c;constraints.ymin]);
end
if isfield(constraints, 'ymin')
    F = [F;C];
    c = abs([c;constraints.ymax]);
end

%% Checking if vector holds NaN values

% Loop through vector and if any of the constraint values are NaN, remove
% that row within the initial matrix and vector

for i = 1:length(c) - 1
    if isnan(c(i))
        F(i,:) = [];
        c(i) = [];
    end
end

%% Loop to find the feasable region and populate the F matrix and c vector with new rows

% Using a while loop and counter to determine the amount of iterations.
% In theory, when all the rows in the following timestep are subsumed by
% the previous set then the iterations can stop as the feasiable region
% cannot get any smaller with the current constraint set. For this, the
% count variable is incremented until its value is equal to that of the
% rows in Fnext. This means that all of the rows have had no effect and can
% therefore stop calculating the next time steps.

F0 = F;
z = length(F0);
count = 0;
i = 1;

while count < z
    % Fnext is the single sided constrant matrix that is generated at the
    % next time step. The values within this matrix are to be evaluated 
    % against the initial single sided constraint matrix row by row. 
    Fnext = F0*(A-B*K)^i;
    i = i + 1;
    count = 0;
    
    for j = 1:z
        
        % Using the linprog function to get the maximum point at which the
        % line hits the feasable region
        x = linprog(-Fnext(j,:),F,c);
        
        % Putting the values of linprog back into the matrix equation to
        % determine if the value is greater than the integer stored in the
        % vector c. 
        const = x(1)*Fnext(j,1,:)+x(2)*(Fnext(j,2,:));
        
        % If the value obtained from const is greater than the value within
        % the vector then the line passes through the feasable region
        % making it smaller. If the value is less than the integer within
        % the vector then it is subsumed by the previous sets.
        if const >= c(j)
            
            % Adding the successful rows from the signle sided contraint
            % matrix and vector into the F matrix and c vector if they
            % are successful
            F = [F;Fnext(j,:)];
            c = [c;c(j)];
        else
            
            % If row does not meet the criteria then a count in incremented
            % If this count reaches the number of rows within Fnext, then
            % it can be said that the feasable region is at its smallest
            % and all preceding lines are subsumed and the process can
            % stop
            count = count + 1;
        end
    end
end

end

