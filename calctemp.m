function [x, t, u] = calctemp(tmax, nt, xmax, nx, method, prop, timeData, tempData, varargin)
% Function for modelling temperature in a space shuttle tile
% D N Johnston  05/01/23 
% Modified by GROUP 44 ME20021  20/04/23 
%
% Input arguments:
% tmax   - maximum time (s)
% nt     - number of timesteps
% xmax   - total thickness (m)
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% timeData - time vector for surface temperatures (s)
% tempData - surface temperature vector (C or K)
%
% Return arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = calctemp(4000, 501, 0.05, 21, 'forward', timeData, tempData);
%

% Set material properties and derived values (LI-900)
% Obtained from NASA document: Structures and Materials: Space Shuttle Tiles, Grades 5-8 - NASA
% Note that we're assuming constant properties.
thermCon = cell2mat(prop(2)); % W/m K; 0.028 BTU/ft/hr/F at 70F and 1 atm
density  = cell2mat(prop(3));    % kg/m^3; 9 lb/ft^3
specHeat = cell2mat(prop(4));    % J/kg/K; 0.15 Btu/lb/F
Bound = prop(1);   % Boundary condition

% Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
u = zeros(nt, nx);
Tinitial = -156 + 275;

alpha = thermCon/(density*specHeat);
p = alpha * dt / dx^2;

% Use interpolation to get outside temperature at time vector t 
% and store it as left-hand boundary vector L.

L = interp1(timeData, tempData, t, "linear", "extrap");

% set initial conditions equal to boundary temperature at t=0.
u(1,:) = Tinitial; % Temperature of space

if strcmp(Bound,'Dirichlet')
    u(:,nx) = Tinitial; % Set inside temperature assuming dirichlet boundary

    % set up index vectors
    i = 2:nx-1;
    im = 1:nx-2; % or i-1
    ip = 3:nx;   % or i+1
elseif strcmp(Bound,'Neumann')
    % set up index vectors assuming Neumann boundary (zero heat flow)
    i = 2:nx;
    im = 1:nx-1; % or i-1
    ip = [3:nx nx-1];   % or i+1
end

u(:, 1) = L; % Outside boundary condition

% Select method and run simulation.
switch method
    case 'Forward differencing'
        % step through time
        for n=1:nt-1
            % calculate internal values using forward differencing
            u(n+1,i) = (1 - 2 * p) * u(n,i) + p * (u(n,im) + u(n,ip));
        end
    case 'DuFort-Frankel'
        % step through time
        for n=1:nt-1
            
            % set index for 'old' point
            if n == 1
                nminus1 = 1; % at first timestep, old point doesn't exist as n-1 = 0.
                             % Use value at timestep 1 instead.
            else
                nminus1 = n-1; % after first timestep, proceed normally.
            end

            % calculate internal values using DuFort-Frankel method
            u(n+1,i) = ((1 - 2*p) * u(nminus1,i) + 2*p * (u(n,im) + u(n,ip)))/(1 + 2*p);
        end
    case 'Backward differencing'
        iVec = 2:nx-1; % set up index vector
        
        % now loop through time
        for n=1:nt-1
            
            % calculate internal values using backward differencing
            b(1)    = 1;    % Could be taken out of loop
            c(1)    = 0;
            d(1)    = L(n); % or u(n,1) (for neumann)
            a(iVec) = -p;
            b(iVec) = 1 + 2*p;
            c(iVec) = -p;
            d(iVec) = u(n,iVec);

            if strcmp(Bound,'Dirichlet')
                a(nx)   = 0;        %
                b(nx)   = 1;        % Inside boundary condition using Dirichlet boundary.
                d(nx)   = u(n,nx);  %
            elseif strcmp(Bound,'Neumann')
                a(nx)   = -2*p;     %
                b(nx)   = 1 + 2*p;  % Inside boundary condition using Neumann boundary.
                d(nx)   = u(n,nx);  %
            end
            
            u(n+1,:) = tdm(a,b,c,d);
        end
    case 'Crank-Nicolson'
        iVec = 2:nx-1; % set up index vector
        
        % now loop through time
        for n=1:nt-1

            % calculate internal values using Crank-Nicolson
            b(1)    = 1;    % Could be taken out if loop
            c(1)    = 0;
            d(1)    = L(n+1); % or u(n+1,1)
            a(iVec) = -p/2;
            b(iVec) = 1 + p;
            c(iVec) = -p/2;
            d(iVec) = (p/2)*u(n,iVec-1) + (1-p)*u(n,iVec) + (p/2)*u(n,iVec+1);

            if strcmp(Bound,'Dirichlet')
                a(nx)   = 0;        %
                b(nx)   = 1;        % Inside boundary condition using Dirichlet boundary.
                d(nx)   = u(n,nx);  %
            elseif strcmp(Bound,'Neumann')
                a(nx)   = -p;       %
                b(nx)   = 1 + p;    % Inside boundary condition using Neumann boundary.
                d(nx)   = p * u(n,nx-1) + (1 - p) * u(n,nx); %
            end
            
            u(n+1,:) = tdm(a,b,c,d);
        end

    otherwise
        error (['Undefined method: ' method])
end

% Produces plots used for troubleshooting.
if nargin > 8
    figure
    % contour plot
    surf(x,t,u)
    % comment out the next line to change the surface appearance
    shading interp  
    
    % Rotate the view
    view(140,30)
    
    %label the axes
    xlabel('\itx\rm - m')
    ylabel('\itt\rm - s')
    zlabel('\itu\rm - deg C')
    title(method)

    if cell2mat(varargin) > 1
        figure
        plot(t,[u(:,1) u(:,nx)], LineWidth=1)
        xlabel('Time (s)')
        ylabel('Temperature')
        grid on
        grid minor
        legend('Outside temperature', 'Inside temperature')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tri-diagonal matrix solution 
function x = tdm(a,b,c,d)
n = length(b);

% Eliminate a terms
for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

x(n) = d(n) / b(n);

% Loop backwards to find other x values by back-substitution
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end