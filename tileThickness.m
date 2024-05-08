function [xTile] = tileThickness(tmax, nt, nx, method, prop, timeData, tempData)
%TILETHICKNESS Determine the optimal tile thickness using shooting method
% 
%   [XTILE] = TILETHICKNESS(TMAX, NT, NX, METHOD, BOUND, TIMEDATA,
%   TEMPDATA) outputs the optimal thickness Xtile

% Target temperature and initial guesses
% Obtained from NASA document: Structures and Materials: Space Shuttle Tiles, Grades 5-8 - NASA
maxTemp = 176 + 273;
guess1 = 0.025;
guess2 = 0.127;

er = maxTemp;
i = 1;
while abs(er) > 1 && i < 100
    [~,~,temp1] = calctemp(tmax, nt, guess1, nx, method, prop, timeData, tempData);
    er1 = max(temp1(:,end)) - maxTemp;
    [~,~,temp2] = calctemp(tmax, nt, guess2, nx, method, prop, timeData, tempData);
    er2 = max(temp2(:,end)) - maxTemp;
    
    % Calculate third guess with shooting method formula
    xTile = guess2 - er2*((guess2 - guess1)/(er2 - er1));
    [~,~,temp] = calctemp(tmax, nt, xTile, nx, method, prop, timeData, tempData);
    er = max(temp(:,end)) - maxTemp;

    guess1 = guess2;
    guess2 = xTile;
    i = i + 1;
end