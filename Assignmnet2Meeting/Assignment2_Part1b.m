%% Assignmnet 2 Part 1 b) - Andrew Paul 100996250
% The second section of this assignment is similar to the first section but
% different boundary conditions are given to make the solution a bit more
% complex. The same finite differences method is used as in part 1 a), but
% now the boundaries of the y axis are set to have a voltage of zero and
% the boudndaries of the x axis have a voltage of V0.

clear

% Set the length and width of the grid and the voltage V0

L = 3;
W =2;
V0 = 1;

% Spacing and number of points of grid

dx = 0.1;
dy = 0.1;
nx = L/dx;
ny = W/dy;

% Parameters of Laplace equation

p1 = 1/(dx^2);
p2 = 1/(dy^2);
p3 = -2*(1/dx^2 + 1/dy^2);

% Generate "G" Matrix

G = zeros(nx*ny,nx*ny);

for i = 2:nx-1
    for j = 2:ny-1
        n = i + (j-1)*nx;
        nym = i + (j-2)*nx;
        nyp = i + j*nx;
        nxm = (i-1) + (j-1)*nx;
        nxp = (i+1) + (j-1)*nx;
        
        G(n,n) = p3;
        G(n,nxm) = p1;
        G(n,nxp) = p1;
        G(n,nym) = p2;
        G(n,nyp) = p2;
    end
end

% Generate "F" Matrix

F = zeros(nx*ny,1);

for i = 1:nx
    n = i;
    G(n,n) = 1;
    
    n1 = i + (ny-1)*nx;
    G(n1,n1) = 1;
end
    
for j = 1:ny
    n = 1 + (j-1)*nx;
    G(n,n) = 1;
    F(n) = V0;
    
    n1 = nx + (j-1)*nx;
    G(n1,n1) = 1;
    F(n1) = 1;  
end

% Solution for F matrix is zero at all corners

F(1) = 0;
F(1 + (ny-1)*nx) = 0;
F(nx) = 0;
F(nx + (ny-1)*nx) = 0;

% Finding solution and reshaping the transpose

A = G\F;
solution = reshape(A,[],ny)';

xrange = linspace(0,L,nx);
yrange = linspace(0,W,ny);

figure(1)
surf(xrange,yrange,solution)
xlabel('x')
ylabel('y')
zlabel('Voltage')
title('Numerical Finite Differences Solution')

analytic = zeros(ny,nx);
steps = 200;

% Create replicas of x and y axes to generate analytic solution

xnew = repmat(linspace(-L/2,L/2,nx),ny,1);
ynew = repmat(linspace(0,W,ny),nx,1)';

% Loop through analytical series solution in steps of 2

for n = 1:2:steps
    analytic = analytic + 1./n.*cosh(n.*pi.*xnew./W)./cosh(n*pi*L/W).*sin(n.*pi.*ynew./W);
end

analytic = analytic*4*V0/pi;

figure(2)
surf(xrange,yrange,analytic)
xlabel('x')
ylabel('y')
zlabel('Voltage')
title('Analytical Finite Differences Solution')

%% 
% After further analyzing the analytical solution it was found that the
% series seems to converage after the first 5 or 6 iterations and shows
% little change as the series continutes.
%
% It is clear that the analytical solution is much easier to impliment as
% it only requires a few lines of code. The situation used for this
% assignmnet is relatively basic and more complicated solutions would
% require a much more complex series solution in order to accomodate
% different boundary conditions.
%
% As for the mesh size, the mesh used in the plots shown are for a dx and
% dy value of 0.1 but when a larger dx and dy value were used it was found
% that the error was much greater. The error was determined by taking an
% average of the difference between the numerical and analytical solutions
% for the total number of itterations.
%
% To conclude, the anayltical solution is a good option for simple cases
% but the numerical solution will be more accurate and potentially easier
% to impliment when facing complicated systems.

        