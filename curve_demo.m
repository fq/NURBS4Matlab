clear all;
% Making java classes available to MATLAB
% Put "nurbs.jar" in current directory
javaaddpath('.\nurbs.jar');
import fq.geom.nurbs.*;

global crv;

%%
% define N control points, cp[N][2]
cp=[32    65     0     1;
    40    52     0     1;
    50    68     0     1;
    62    69     0     1;
    63    79     0     1;
    52    84     0     1;
    61    99     0     1;
    70   100     0     1;
    73    88     0     1;
    87    63     0     1;
    94    64     0     1;
   106    68     0     1];

degree = 6;

% create a B-spline curve
crv = NurbsCurve(cp,degree);

% USAGE:
% Function to calculate a point(or points) of this curve at a parameter value
%   pt = crv.PointAt(u)
%   x = pt(1); y = pt(2); z = pt(3);
%
% Function to calculate derivates of this curve at a parameter value
%   ders = crv.DerivsAt(u,2);
%   x = ders(1,1);y = ders(1,2);
%   dx = ders(2,1);dy = ders(2,2);
%   dx2 = ders(3,1);dy2 = ders(3,2);

NN = 100;
% sample code
uu = linspace(0,1,NN); % generate a set of parameter values
pt = crv.PointAt(uu); % calculate points
tt = crv.TangentAt(uu);  % calculate unit tangent vectors
nt = crv.NormalAt(uu); % calculate unit proncipal normal vectors
k = crv.CurvatureAt(uu); % calculate curvature
r = 1./k;

figure;
hold off;
%plot3(cp(:,1),cp(:,2),cp(:,3),'-ro');
plot(cp(:,1),cp(:,2),'-ro');
hold on;
%plot3(pt(:,1),pt(:,2),pt(:,3),'-','LineWidth',2);
plot(pt(:,1),pt(:,2),'-','LineWidth',2);
%quiver3(pt(:,1),pt(:,2),pt(:,3),nt(:,1),nt(:,2),nt(:,3));
quiver(pt(:,1),pt(:,2),nt(:,1),nt(:,2));

title('NURBS curve with unit principal normal vectors');
%axis equal vis3d;
