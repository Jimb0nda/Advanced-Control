function [Aa, Ba, Ca] = matrix_aug(A,B,C)
% MATRIX AUGMENTATION 
% Create augmented matrices for MPC computation
%
% Aa = [A B; 0 I]    - where 0 is rxn and I is rxr
% Ba = [B; I]        - where I is rxr
% Ca = [C 0]         - where 0 is mxr

%% Augment Matrices

Aa = [A B; zeros(size(B,2), size(A,1)) eye(size(B,2))];

Ba = [B; eye(size(B,2))];

Ca = [C zeros(size(C,1), size(B,2))];

end

