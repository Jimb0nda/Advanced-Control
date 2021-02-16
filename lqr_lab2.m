%% Matrices set up

A = [0 1 0 0; 20.601 0 0 0; 0 0 0 1; -0.4905 0 0 0];
B = [0;-1;0;0.5];
Q = [100 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R = 1;


%% System response to the following initial conditions

x0 = [0.1;0;0;0];

[K,P] = lqr(A,B,Q,R);

figure(1);
initial(feedback(ss(A,B,eye(4),[]),K),x0,10)
damp(feedback(ss(A,B,eye(4),[]),K))

%% Testing observability to determine which state is observable

C1 = [1 0 0 0];
C2 = [0 1 0 0];
C3 = [0 0 1 0];
C4 = [0 0 0 1];

obvs1 = [C1; C1*A; C1*A^2; C1*A^3];
obvs2 = [C2; C2*A; C2*A^2; C2*A^3];
obvs3 = [C3; C3*A; C3*A^2; C3*A^3];
obvs4 = [C4; C4*A; C4*A^2; C4*A^3];

s_ob1 = rank(obvs1)
s_ob2 = rank(obvs2)
s_ob3 = rank(obvs3)
s_ob4 = rank(obvs4)