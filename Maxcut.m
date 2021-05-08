clc; clear;
% Initial conditions
n = 5; %n x n set of points (n=30 is 900 points)
J = -1;

%Sets up the adjacency matrix for the 2d ising model
A = zeros(n^2);
for i = 1:n^2-1
    if (mod(i,n) ~= 0)
        A(i,i+1)=-J;
        A(i+1,i)=-J;
    end
end
for i = 1:n^2-n
    A(i,i+n)=-J;
    A(i+n,i)=-J;
end

%Calculates the maximum cut and displays the cut
[cut,y] = gw_MaxCut(A, 1000);
fprintf("Value of found max-cut = %d\n",cut)
fprintf("Value of energy = %d\n", -1*cut/2)
y=reshape(y,n,n);
hold on
title("Goemans-Williamson",'FontSize',18)
imshow(y,'InitialMagnification',1500) %May have to change magnification size for larger n

function [cut, y] = gw_MaxCut(A, T)
    [n,~] = size(A);
    % Use CVX
    cvx_begin quiet
        variable X(n,n) symmetric
        minimize trace(A*X)
            diag(X) == ones(n,1);
            X == semidefinite(n);
    cvx_end
    
    % Factor the matrix
    U = chol(X);
    
    % Initialize the values
    cut = zeros(T,1);
    y = zeros(n,T);
    
    % Perform T trials
    for i = 1:T
        % Randomly select a hyperplane
        r = mvnrnd(zeros(n,1),diag(ones(n,1)))';
        % Select which group each node falls into
        y(:,i) = sign(U*r);
        % Calculate the value of the cut
        cut(i) = (sum(A(:)) - y(:,i)'*A*y(:,i))/4;
    end
    % Take the largest of the trials as our answer.
    [cut,i]=max(cut);
    y = y(:,i);
end