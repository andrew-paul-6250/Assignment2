%% Assignment 2 Part 2 - Andrew Paul 100996250
% The final section of this assignment also involeved using the finite
% difference method but to model the current flow of a rectangular region
% with limiting boxes. The boxes were given some conductivity value and the
% effect of the "bottle-neck" effect on the current was investigated.

clear

nx = 50;
ny = 50;

% Create sparse G matrix
G = sparse(nx*ny,nx*ny);

% Conductivity outside box
sigma1 = 1;
% Conductivity inside box
sigma2 = 10^-2;

% Generate F matrix to set boundary conditions
F = zeros(nx*ny,1);

% Change for difference in bottle neck width
Lb = 0.4;
Wb = 0.6;

% Create matrix for mapping the conductivity and loop through to assign
% conductivity values for the given conditions
condMap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        if (i>=Lb*nx && i<=Wb*nx && j<=Lb*ny) || (i>=Lb*nx && i<=Wb*nx && j>=Wb*ny)
            condMap(i,j) = sigma2;
        else
            condMap(i,j) = sigma1;
        end
    end
end

% Loop through to set boundary conditions

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
            
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            
        elseif j == 1
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nyp = j+1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            
        elseif j == ny
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
        
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            
        else 
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nym = j-1+(i-1)*ny;
            nyp = j+1+(i-1)*ny;
            
            rxm = (condMap(i,j) + condMap(i-1,j))/2;
            rxp = (condMap(i,j) + condMap(i+1,j))/2;
            rym = (condMap(i,j) + condMap(i,j-1))/2;
            ryp = (condMap(i,j) + condMap(i,j+1))/2;
        
            G(n,n) = -(rxm+rxp+ryp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
        
    end
    
end

% Find voltage values using matrix operations
V = G\F;

% Create matrix to map voltage and loop through matrix to assign values
% from calculated voltage matrix
VMap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        VMap(i,j) = V(n);
    end
end

% Plots

% Conductivity plot
figure(1)
surf(condMap)
title('Conductivity of region')
xlabel('y')
ylabel('x')
zlabel('Conductivity (ohm meters)')
view(0,90)

%% 
% The plot above shows the the conductivity has one constant value inside the
% boxes and a constant value outisde of the boxes which is expected as
% these were the conditions set to the system.


% Voltage plot
figure(2)
surf(VMap)
title('Voltage distribution over region')
xlabel('y')
ylabel('x')
colorbar
view(0,90)

%% 
% The voltage plot shown above gives the expected distribution as one side
% of the region has a voltage of V0 where the other side is set to zero.
% The barrier in the middle blocks uniform distribution to the other side
% of the region which is expected. As the conductivity inside the boxes is
% less than the conductivity outside the barrier boxes.

% Gradient used to plot electric field lines
[Ex,Ey] = gradient(VMap);

% Electric field plot
figure(3)
quiver(Ex,Ey)
title ('Electric Field Lines')
xlabel('x')
ylabel('y')
xlim([0 50])
ylim([0 50])

%% 
% The plot above shows that the electric field lines are strongest between
% the areas which have a lower conductivity. This is execpted as it is
% similar to the model of a parallel plate capcitor which has two plates of
% a larger conductivity seperated by a region with a lower conductivity
% creating a stronger electric field between the two higher conductivity
% regions.


% Calculation of current density
jx = condMap.*Ex;
jy = condMap.*Ey;

% Current density plot
figure(4)
quiver(jx, jy)
title('Current Density Inenstiy Lines')
xlim([0 50])
ylim([0 50])

%% 
% The current density plot shows the largest current density is between the
% boxes which is expected as there is the same amount of current being
% squeezed through a smaller region.


% Calculation of current flow
Ex = -Ex;
Ey = -Ey;

xflow = condMap.*Ex;
yflow = condMap.*Ey;

x0sum = sum(xflow(:,1));
xsum = sum(xflow(:,nx));
y0sum = sum(yflow(:,1));
ysum = sum(yflow(:,nx));

xCurrent = (x0sum + xsum)/2;
yCurrent = (y0sum + ysum)/2;
% output the current for each condition
totCurrent = sqrt(xCurrent^2 + yCurrent^2)

meshMultiple = [1 2 3 4 5];
% mesh current was cacluclated using different mesh multiples
meshCurrent = [0.1628 0.1719 0.175 0.1765 0.1775];

figure(5)
plot(meshMultiple, meshCurrent)
title('Current vs Mesh Densities')
xlabel('Mesh density')
ylabel('Current (A)')

%% 
% The plot above shows how the current saturates as the mesh becomes finer
% with a multiplier, this is expected as a finer mesh grid will allow for a
% more accurate solution of how the current is responding which will
% ultimaltly converge to be one finite value.


width = [0.4 0.2 0.16 0.12 0.08 0.04];
% different current values were calculated using different bottle-neck
% widths
widthCurrent = [0.255 0.1628 0.1442 0.1344 0.101 0.0707];

figure(6)
plot(width, widthCurrent, '-o')
title('Current vs Bottle-neck Width')
xlabel('Bottle-neck Width')
ylabel('Current (A)')

%% 
% The figure above displays what happens to the current as the bottle neck
% becomes smaller, as expected the current gets smaller as the bottle neck
% becomes smaller as it is restricting the amount of current which can flow
% through the region. More points could be added to show a better trend
% which would begin to satruate.


cond = [0.2 0.5 1 2 3];
% different current values were calculated using different conductivities
condCurrent = [0.0576 0.0985 0.1628 0.2897 0.4162];

figure(7)
plot(cond, condCurrent, '-o')
title('Current vs Conductivity')
xlabel('Conductivity (ohm meter)')
ylabel('Current (A)')

%% 
% Finally, the figure above shows that the conductivity is inversely
% proportional to the current and thus the current decreases when the
% conductivity increases and increases when the conductivity decreases.
% This is expected 

            


