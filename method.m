function [ v,pres,rho ] = method( X,T,XP,TP,referenceElement,mu,nu,confined,nStep,dom,dt,nx,ny,periodic )

% Matrices arising from the discretization
% MV is the mass matrix for velocity field
% Mh is the mass matrix for rho


[MV,Mh,Kh,G,~] = FEM_system(X,T,XP,TP,referenceElement);


[ndofP,ndofV] = size(G); 

if confined
   nunkP = ndofP-1;
   disp(' ')
   disp('Confined flow. Pressure on lower left corner is set to zero');
   G(1,:) = [];
else
   nunkP = ndofP;
end

pres = zeros(nunkP,nStep+1);
v=zeros(size(X,1)*2,nStep+1);
rho=zeros(size(X,1),nStep+1);

% set initial velocity and density field
[ v_0,rho_0 ] = InitialCondition( X );

%velocity
v(:,1)=v_0;

%density
rho(:,1) = rho_0;

%constant langrange multipliers for rho

[Ah,bh,nDirh]=BC_h(X,dom);


% Matrix and r.h.s vector to impose Dirichlet boundary conditions using
% Lagrange multipliers
for n=1:nStep
 
 %update source term for density field in each time step
 [Ch ,fh]= ConvSource(X,T,referenceElement,v(:,n));
 
  Arho = Mh/dt + 0.5*(Ch+Kh*mu);
  brho = -(Ch+Kh);
 
 Ktoth=[Arho Ah';Ah zeros(nDirh,nDirh)];
 
 [Lh,Uh] = lu(Ktoth);
 
 ftoth=[(brho*rho(:,n) + fh);bh];
 
 solh=Uh\(Lh\ftoth);
 Du=solh(1:length(fh));
 rho(:,n+1)=rho(:,n)+Du;
 
 % Until now, we have advanced the rho in one time step,now we could
 % advance velocity field
 
 % current time
 tn=n*dt;
 % non-constant langrangemultiplier for V, updated in each time step
[A_DirV, b_DirV, nDirV] = BC_V(X,dom,ndofV,tn,periodic);

%update stiffness matrix for Velocity field in each time step
 K = Stiff_nu( X,T,referenceElement,rho(:,n+1),nu );

 %  Chorin-Temam Method

%first step 
% lagrange multiplier to apply Dirichlet boudary condition
Ktot1=[MV+dt*K A_DirV';A_DirV zeros(nDirV)];
btot1=[MV*v(:,n);b_DirV];
Ktot1=sparse(Ktot1);
btot1=sparse(btot1);

%lu factorization
[LV_int, UV_int]=lu(Ktot1);

sol_int=UV_int\(LV_int\btot1);

v_int=sol_int(1:size(X,1)*2);

% second step
% viscosity splitting fractional-step method Blasco, Codina and Huerta

 Ktot2 = [(MV/dt)+K        A_DirV'             G'
            A_DirV    zeros(nDirV,nDirV)     zeros(nDirV,nunkP)
            G          zeros(nunkP,nDirV)    zeros(nunkP,nunkP)];
 btot2=[MV*v_int/dt+K*v_int;b_DirV;zeros(nunkP,1)];
 Ktot2=sparse(Ktot2);
 btot2=sparse(btot2);
 [LV2,UV2]=lu(Ktot2);
 solv=UV2\(LV2\btot2);
 
 v(:,n+1)=solv(1:size(X,1)*2);
 pres(:,n+1)=solv(size(X,1)*2+nDirV+1:size(X,1)*2+nDirV+nunkP);
 


end


if confined
   pres=[zeros(1,nStep+1);pres]; 
end


end













