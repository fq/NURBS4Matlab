clear all;
% making java classes available to MATLAB
% Put "nurbs.jar" in current directory
javaaddpath('.\nurbs.jar');
import fq.geom.nurbs.*;

%%
% define U x V control points, cp[U][V][4]
cp(:,:,1) = [ 0 0 0 0 0 ;
              1 1 1 1 1 ;
              2 2 2 2 2 ;
              3 3 3 3 3 ;
              4 4 4 4 4]; % x coordinates
cp(:,:,2) = [ 0 1 2 3 4 ;
              0 1 2 3 4 ;
              0 1 2 3 4 ;
              0 1 2 3 4 ;
              0 1 2 3 4 ]; % y coordinates
cp(:,:,3) = [ 0 -1 0 -1 0;
              -1 -2 -1 -2 -1;
              0 0 0 0 0 ;
              1 2 1 2 1 ;
              0 1 0 1 0]; % z coordinates
cp(:,:,4) = [ 1 1 1 1 1 ;
              1 1 1 1 1 ;
              1 1 1 1 1 ;
              1 1 1 1 1 ;
              1 1 1 1 1 ]; % weights

degreeU = 3;degreeV=3;
srf = NurbsSurface(cp,degreeU,degreeV);

% sample code to draw the surface
UU =31;VV=31;
u = linspace(0,1,UU);
v = linspace(0,1,VV);
p = srf.PointAt(u,v);
x = p(:,:,1);y=p(:,:,2);z=p(:,:,3);
n = srf.NormalAt(u,v);
nx = n(:,:,1);ny=n(:,:,2);nz=n(:,:,3);

surf(x,y,z);
axis equal vis3d
hold on;
quiver3(x,y,z,nx,ny,nz);
title('NURBS surface with unit normal vectors');
