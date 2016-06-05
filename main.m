% This program solves Stokes problem in a square domain


clear; close all; clc

addpath('Func_ReferenceElement')

dom = [0,2,0,3]; 


% Element type and interpolation degree
% (0: quadrilaterals, 1: triangles, 11: triangles with bubble function)
%elemV = 0; degreeV = 2; degreeP = 1;
 elemV = 1; degreeV = 2; degreeP = 1;
% elemV = 11; degreeV = 1;  degreeP = 1; 
if elemV == 11
    elemP = 1; 
else
    elemP = elemV; 
end
referenceElement = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP); 

hh = cinput('Spatial mesh size',0.2);
nx=(dom(2)-dom(1))/hh;
ny=(dom(4)-dom(3))/hh;
adapted = 0;

[X,T,XP,TP] = CreateMeshes(dom,nx,ny,referenceElement,adapted);

figure; PlotMesh(T,X,elemV,'b-');
figure; PlotMesh(TP,XP,elemP,'r-');

% Time discretization
disp('');
tEnd = cinput('End time', 2*pi);
disp('');
nStep = cinput('Numnber of time-steps', 110);
dt = tEnd / nStep; 
%Courant = sqrt(Conv(1)^2+Conv(2)^2)*dt/hh; 
%disp(['Courant number: ', num2str(Courant)]); 


disp('Confined flow:');
disp('   [1] Yes');
disp('   [2] No');
confined=cinput('Pressure on lower left corner is set to zeros',1);
disp('');
nu=cinput('Diffusion coefficient nu for velocity',0.02);
disp('');
mu=cinput('Diffusion coefficient mu for rho',1);
disp('');
periodic=cinput('Periodic time for Dricihlet boundary in velocity field',2);
[ v,pres,rho ] = method( X,T,XP,TP,referenceElement,mu,nu,confined,nStep,dom,dt,nx,ny,periodic );





%% postprocess


xx  = reshape(X(:,1), degreeV*nx+1, degreeV*ny+1)'; 
yy  = reshape(X(:,2), degreeV*nx+1, degreeV*ny+1)';

%figure(3)

density=moviein(size(rho,2));

for n=1:size(rho,2)

figure(3)
clf(3,'reset');
rhoo = reshape(rho(:,n), degreeV*nx+1, degreeV*ny+1)';  
surface(xx,yy,rhoo,'FaceColor','interp');
set(gca,'FontSize',16)
grid on
view(3)
density(:,n)=getframe(gcf);
pause(0.1)

end

%movie2avi(density,'density','compression','none');

pressure=moviein(size(pres,2));
for n=1:size(pres,2)
    
 if degreeP == 0
     PlotResults(X,T,pres(:,n),referenceElement.elemP,referenceElement.degreeP)
 else
    PlotResults(XP,TP,pres(:,n),referenceElement.elemP,referenceElement.degreeP)
 end
 pressure(:,n)=getframe(gcf);
 pause(0.1)
end

%movie2avi(pressure,'pressure','compression','none');


velocityfield=moviein(size(v,2));
for n=1:size(v,2)

velo=[v(1:2:end-1,n),v(2:2:end,n)];
nPt = size(X,1); 
figure(5);
clf(5,'reset');
quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1),velo(1:nPt,2));
hold on 
plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
axis equal; axis tight

velocityfield(:,n)=getframe(gcf);
 pause(0.1)
end
%movie2avi(velocityfield,'velocityfield','compression','none');

streamline=moviein(size(v,2));
for n=1:size(v,2)

    
    
velo=[v(1:2:end-1,n),v(2:2:end,n)]; 
    
PlotStreamlines(X,velo,dom);
streamline(:,n)=getframe(gcf);

pause(0.1)

end
%movie2avi(streamline,'streamline','compression','none');