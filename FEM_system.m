function [MV,Mh,Kh,G,f] = FEM_system(X,T,XP,TP,referenceElement)
% [K,G,f] = Stokes_system(X,T,XP,TP,referenceElement)
% Matrices K, G and r.h.s vector f obtained after discretizing a Stokes problem
%
% X,T: nodal coordinates and connectivities for velocity
% XP,TP: nodal coordinates and connectivities for pressure
% referenceElement: reference element properties (quadrature, shape functions...)


elem = referenceElement.elemV;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
NP = referenceElement.NP; 
ngeom = referenceElement.ngeom; 

% Number of elements and number of nodes in each element
[nElem,nenV] = size(T);
nenP = size(TP,2); 

% Number of nodes
nPt_V = size(X,1);
if elem == 11
    nPt_V = nPt_V + nElem; 
end
nPt_P = size(XP,1);

% Number of degrees of freedom 
nedofV = 2*nenV; 
nedofP = nenP;
ndofV = 2*nPt_V; 
ndofP = nPt_P; 


MV= zeros(ndofV,ndofV);
Kh= zeros(nPt_V,nPt_V);
Mh= zeros(nPt_V,nPt_V);
G = zeros(ndofP,ndofV); 
f = zeros(ndofV,1);


% Loop on elements
for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
    TPe = TP(ielem,:); 
    % Coordinates of the nodes in element ielem
    Xe = X(Te(1:ngeom),:);
    % Degrees of freedom in element ielem
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofV);
    TPe_dof = TPe; 
    
    % Element matrices
    [MVe,Mhe,Khe,Ge,fe] = EleMatStokes(Xe,ngeom,nedofV,nenV,nedofP,ngaus,wgp,N,Nxi,Neta,NP);
    
    % Assemble the element matrices
  
    MV(Te_dof, Te_dof) = MV(Te_dof, Te_dof) + MVe;
    Kh(Te,Te)=Kh(Te,Te)+Khe;
    Mh(Te,Te)=Mh(Te,Te)+Mhe;
    G(TPe_dof,Te_dof) = G(TPe_dof,Te_dof) + Ge; 
    f(Te_dof) = f(Te_dof) + fe;
   
end






function [MVe,Mhe,Khe,Ge,fe] = EleMatStokes(Xe,ngeom,nedofV,nenV,nedofP,ngaus,wgp,N,Nxi,Neta,NP)
%

Khe= zeros(nenV,nenV);
MVe = zeros(nedofV,nedofV);
Mhe= zeros(nenV,nenV);
Ge = zeros(nedofP,nedofV);
fe = zeros(nedofV,1);



% Loop on Gauss points 
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    NP_ig = NP(ig,:); 
    Jacob = [
        Nxi_ig(1:ngeom)*(Xe(:,1))	Nxi_ig(1:ngeom)*(Xe(:,2))
        Neta_ig(1:ngeom)*(Xe(:,1))	Neta_ig(1:ngeom)*(Xe(:,2))
        ];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    nx = res(1,:);
    ny = res(2,:);
    
	Ngp = [reshape([1;0]*N_ig,1,nedofV); reshape([0;1]*N_ig,1,nedofV)];
    % Gradient
    Nx = [reshape([1;0]*nx,1,nedofV); reshape([0;1]*nx,1,nedofV)];
    Ny = [reshape([1;0]*ny,1,nedofV); reshape([0;1]*ny,1,nedofV)];
    % Divergence
    dN = reshape(res,1,nedofV);

  
    Khe = Khe+(nx'*nx+ny'*ny)*dvolu;
    MVe=MVe+Ngp'*Ngp*dvolu;
    Mhe=Mhe+N_ig'*N_ig*dvolu;
    Ge = Ge - NP_ig'*dN*dvolu;
    
    x_ig = N_ig(1:ngeom)*Xe; 
    f_igaus = SourceTerm1(x_ig);
    fe = fe + Ngp'*f_igaus*dvolu; 
end

