function [Err] = Analyser(tmax, nt, xmax, nx, method, prop, unit, timeData, tempData)
%ANALYSER Analysis of time and spatial step.
% 
%   [ERR] = ANALYSER(TMAX, NT, XMAX, NX, METHOD, BOUND, UNIT, TIMEDATA, 
%   TEMPDATA) outputs an estimation of the percentage error due to number 
%   of step NX. The temperature at time TMAX is plotted against timestep 
%   for each method and the temperature at a timestep of 0s is linearly
%   extrpolated to give a estimation of the analytical solution. The 
%   percentage error of temperature due to step size nx is then calculated 
%   and output.
%   The inputs TMAX, NT, XMAX, NX, METHOD, BOUND, UNIT, TIMEDATA, TEMPDATA 
%   are the same as those used in the function calctemp().

% Find the last temperature at the inner boundary
[~, ~, u] = calctemp(tmax, nt, xmax, nx, method, prop, timeData, tempData);
temp = u(end, nx);

% Loop through a large range of timesteps for each method and store the
% final temperature of the inside boundary
j = 1;
for nt = 701:20:1501
    dt(j) = tmax/(nt-1);
    % disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s'])
    [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Forward differencing', prop, timeData, tempData);
    uf(j) = u(end, nx);
    [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Backward differencing', prop, timeData, tempData);
    ub(j) = u(end, nx);
    [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'DuFort-Frankel', prop, timeData, tempData);
    ud(j) = u(end, nx);
    [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Crank-Nicolson', prop, timeData, tempData);
    uc(j) = u(end, nx);
    j = j + 1;
end

% Calculate the gradient of each method.
gradF = polyfit(dt,uf,1);
gradB = polyfit(dt,ub,1);
gradD = polyfit(dt,ud,1);
gradC = polyfit(dt,uc,1);

% Average teperature at 0s timestep
sol = mean([gradF(2) gradB(2) gradD(2) gradC(2)]);

% Calculate error due to timestep
switch unit
    case "K"
        Err = abs(temp - sol) * (100/(sol));
    case "°C"
        Err = abs(temp - sol) * (100/(sol - 273.15));
    case "°F"
        Err = abs(temp - sol) * (100/((sol - 273.15) * 9/5 + 32));
end

% % The same as above but with spacial step
% j = 1;
% for nx = 3:1:81
%     dx(j) = xmax/(nx-1);
%     % disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s'])
%     [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Forward differencing', Bound, timeData, tempData);
%     uff(j) = u(end, nx);
%     [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Backward differencing', Bound, timeData, tempData);
%     ubb(j) = u(end, nx);
%     [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'DuFort-Frankel', Bound, timeData, tempData);
%     udd(j) = u(end, nx);
%     [~, ~, u] = calctemp(tmax, nt, xmax, nx, 'Crank-Nicolson', Bound, timeData, tempData);
%     ucc(j) = u(end, nx);
%     j = j + 1;
% end

%% PLot the results
% figure
% subplot(2,1,1)
% plot(dt, [uf; ub; ud; uc])
% 
% ylim([380 410])
% legend ('Forward', 'Backward', 'DuFort-Frankel', 'Crank-Nicolson', Location='best')
% ylabel('Temperature (K)')
% xlabel('Timestep (s)')
% grid on
% grid minor
% 
% subplot(2,1,2)
% plot(3:81, [uff; ubb; udd; ucc])
% ylim([385 387])
% legend ('Forward', 'Backward', 'DuFort-Frankel', 'Crank-Nicolson', Location='southeast')
% ylabel('Temperature (K)')
% xlabel('Number of spacial steps')
% grid on
% grid minor