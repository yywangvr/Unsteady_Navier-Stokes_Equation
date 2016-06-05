function [A, b, nDir] = BC_V(X,dom,n,tn,periodic)
% [A, b, nDir, confined] = BC(X,dom,n)
% Matrices to impose Dirichlet boundary conditions using Lagrange
% multipliers on a rectangular domain
% Input: 
%    X: nodal coordinates
%    dom: domain description [x1,x2,y1,y2]
%    n: number of velocity degrees of freedom
% Output:
%    A,b: matrix and r.h.s. vector to impose the boundary conditions using
%         Lagrange multipliers
%    nDir: number of prescribed degrees of freedom
%    confined: 

x1 = dom(1); 
tol = 1e-6; 
 

nodesX1 = find(abs(X(:,1)-x1) < tol); 

 
Vx=DirBC(tn,periodic);

C = [2*nodesX1-1, ones(size(nodesX1))*Vx];

nDir = length(C);  

A = zeros(nDir,n); 
A(:,C(:,1)) = eye(nDir); 

b =C(:,2) ;

      

end


function Vx=DirBC(t,periodic)


T=periodic;

w=2*pi/T;
Vx=1+sin(w*t-pi/2);


end


