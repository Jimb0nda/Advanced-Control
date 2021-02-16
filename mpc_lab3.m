%% Set up

clear

% State Space System
A = [1.9 -0.88; 1 0];
B = [2;0];
C = [-0.5 0.6];

% Prediction Horizon and Control Horizon
Np = 20;
Nc = 5;

% Move Penalty Weight
r = eye(Nc);

% Sample period
T = 0.2;

% Initial Conditions
xa = [0;0;0];

% If desired trajectory S is just a fixed set-point s, then, 
% S = [1 1 ...1^Np-1]^T*s  %s = 1
S = ones(Np, 1);

% Time steps
Ts = 50;

% matrix of time steps to Sample Period T by Ts steps 
t = (1:Ts)*T;

%% Augmented Matrices

Aa = [A B; 0 0 1];
Ba = [B;1];
Ca = [C 0];

%% Contructing F and Phi matrices

F = zeros(Np,3);
for i = 1:Np
    
    F(i,:) = Ca*Aa^i;
    
end

Phi = zeros(Np,Nc);
%Construct first column of Phi matrix
for l = 1:Np
   
    Phi(l,1) = Ca*Aa^(l-1)*Ba;
    
end

%Populate the rest of the Phi matrix
for j=2:Nc
   
    Phi(:,j) = [0;Phi(1:end-1,j-1)];
    
end

%% Solving for the equation du_opt

% du_opt equation, du_opt = (Phi^T*Phi+R)^-1*Phi^T*(S-F*xa[k]);

for i = 1:Ts

    du_opt = inv(Phi.'*Phi+r)*Phi.'*(S-F*xa);
    du(i) = du_opt(1);
    xa = Aa*xa+Ba*du_opt(1);
    u = sum(du);
    u_sum(i) = u; 
    y(i) = Ca*xa;

end

%stairs(t,y)
%hold on
%stairs(t,u_sum)

%% The MPC Toolbox and Terminal Constraints

% Create discrete state space model
sys = ss(A,B,C,[],T);

% Create an MPC object
mpcobj = mpc(sys);

% Change defaults for the mpcobj
mpcobj.PredictionHorizon = Np;
mpcobj.ControlHorizon = Nc;
mpcobj.Weights.ManipulatedVariablesRate=1;

%sim(mpcobj, Ts, 1)

%% Terminal Constraints

Cnew = [-0.5 0.6; 1 -1];

% Create new state space model with new C matrix
sysTC = ss(A,B,Cnew,[],T);

% Create a new MPC object
mpcobj2 = mpc(sysTC);

% Change defaults for the mpcobj
mpcobj2.PredictionHorizon = 5;
mpcobj2.ControlHorizon = 2;
mpcobj2.Weights.ManipulatedVariablesRate=1;
setterminal(mpcobj2, struct('Min', [1 0], 'Max', [1 0]));

sim(mpcobj2, Ts, [1 0])
