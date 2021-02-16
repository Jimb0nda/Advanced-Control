clear

A = [5 4;1 2];
B = [1;0];
C = eye(2);
Q = eye(size(A));
R =[1];

sys = ss(A,B,C,[]);

constr = struct('umin',1,'umax',1,'xmin',[-2; -2],'xmax',[2; 2]);
