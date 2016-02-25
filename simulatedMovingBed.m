function objective = simulatedMovingBed(varargin)

% =============================================================================
% This is the main function which is charge of switching to reach the
% cyclic steady state. The layout of the columns and the configuration is
% listed as follow:
%
%                                       FOUR-ZONE
%              4-column SMB                                       8-column SMB
% Extract                          Feed       |    Extract                           Feed
%       \                          /          |         \                            /
%        --------Zone II(b)--------           |          --------Zone II(c/d)--------
%        |                        |           |          |                          | 
% Zone I(a)                  Zone III(c)      |     Zone I(a/b)               Zone III(e/f)
%        |                        |           |          |                          | 
%        --------Zone IV(d)--------           |          --------Zone IV(h/g)--------
%       /                          \          |         /                            \
% Desorbent                       Raffinate   |   Desorbent                         Raffinate
%
%             12-column SMB                                       16-column SMB
% Extract                            Feed       |    Extract                         Feed
%       \                            /          |         \                          /
%        -- ----Zone II(d/e/f)-------           |          -----Zone II(e/f/g/h)-----
%        |                          |           |          |                        | 
% Zone I(c/b/a)                Zone III(g/h/i)  |  Zone I(a/b/c/d)           Zone III(i/j/k/l)
%        |                          |           |          |                        | 
%        -------Zone IV(l/k/j)-------           |          -----Zone IV(p/o/n/m)-----
%       /                            \          |         /                          \
% Desorbent                         Raffinate   |   Desorbent                       Raffinate
%
%                                       FIVE-ZONE
%              5-column SMB                                       10-column SMB
%    Ext2                          Feed       |      Ext2                            Feed
%       \                          /          |         \                            /
%        --------Zone II(c)--------           |          --------Zone III(e/f)--------
%        |                        |           |          |                           | 
% Zone II(b)                      |           |     Zone II(d/c)                     |
%        |                        |           |          |                           |
% Ext1 --                    Zone IV(d)       |   Ext1 --                        Zone IV(g/h)
%        |                        |           |          |                           |
% Zone I(a)                       |           |     Zone I(b/a)                      |
%        |                        |           |          |                           | 
%        --------Zone V(e)---------           |          ---------Zone V(j/i)---------
%       /                          \          |         /                            \
% Desorbent                       Raffinate   |   Desorbent                         Raffinate
%
%             15-column SMB                                       20-column SMB
%    Ext2                            Feed       |      Ext2                              Feed
%       \                            /          |         \                              /
%        -------Zone II(g/h/i)-------           |          -------Zone III(i/g/k/l)-------
%        |                          |           |          |                             | 
% Zone II(f/e/d)                    |           | Zone II(h/g/f/e)                       |
%        |                          |           |          |                             |
% Ext1 --                    Zone IV(j/k/l)     |   Ext1 --                        Zone IV(m/n/o/p)
%        |                          |           |          |                             |
% Zone I(c/b/a)                     |           | Zone I(d/c/b/a)                        |
%        |                          |           |          |                             | 
%        -------Zone V(o/n/m)--------           |          -------Zone V(t/s/r/q)---------
%       /                            \          |         /                              \
% Desorbent                         Raffinate   |   Desorbent                           Raffinate
%
% Fluid phase goes from Zone I to Zone II to Zone III, while the ports switch direction
% is from Zone I to Zone IV to Zone III;
% =============================================================================


    global string stringSet;

    tTotal = tic;

    [opt, interstVelocity, Feed] = getParameters(varargin{:});

%   Initialize the starting points, currentData
    currentData    = cell(1, opt.nColumn);
    transitionData = cell(1, opt.nColumn);
    tempData       = cell(1, opt.nColumn);
    columnNumber   = cell(1, opt.nColumn);

    for k = 1:opt.nColumn
%       currentData stores the outlets of each interval of columns
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents);
        currentData{k}.lastState = []; 

%       transitionData stores the composition profile of nInterval pieces
        transitionData{k}.outlet.time = linspace(0, opt.switch*opt.nInterval, opt.timePoints*opt.nInterval);
        transitionData{k}.outlet.concentration = zeros(opt.timePoints*opt.nInterval, opt.nComponents); 

%       This is the temporary data for dynamic trajectory plotting
        tempData{k}.concentration = cell(1, opt.nInterval);
    end

%   Number the columns for the sake of plotting
    for k = 1:opt.nColumn
        if k == 1
            columnNumber{1} = opt.nColumn;
        else
            columnNumber{k} = k-1;
        end
    end

% 	Construct the string in order to tell simulator the calculation sequence
	stringSet = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm'...
                 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'...
                 'aa' 'bb' 'cc' 'dd' 'ee' 'ff' 'gg' 'hh' 'ii' 'jj' 'kk' 'll' 'mm'...
                 'nn' 'oo' 'pp' 'qq' 'rr' 'ss' 'tt' 'uu' 'vv' 'ww' 'xx' 'yy' 'zz'};

    if opt.nColumn > length(stringSet)
        error('The simulation of %3g-column case in %3g-zone is not finished so far', opt.nColumn, opt.nZone);
    end

    sequence = cell2struct( columnNumber, stringSet(1:opt.nColumn), 2 );

    string = char(stringSet(1:opt.nColumn));
 
% 	Specify the column for the convergence checking
    if opt.nZone == 4
    	convergIndx = sum(opt.structID(1:2));
	elseif opt.nZone == 5
		convergIndx = sum(opt.structID(1:3));
    end

%   Preallocation
%   The dimension of plotData (columnNumber x switches)
%          t_s   2*t_s  3*t_s  4*t_s
%          {1}    {1}    {1}    {1}
%          {2}    {2}    {2}    {2}
%          {3}    {3}    {3}    {3}
%          {4}    {4}    {4}    {4}
    plotData = cell(opt.nColumn,opt.nColumn);

    if opt.nZone == 4
        dyncData = cell(2, opt.nMaxIter);
    elseif opt.nZone == 5
        dyncData = cell(3, opt.nMaxIter);
    end

%   convergPrevious is used for stopping criterion
    convergPrevious = transitionData{convergIndx}.outlet.concentration;


%-----------------------------------------------------------------------------------------
%   Main loop
    for i = 1:opt.nMaxIter

%		Switching the ports, in the countercurrent manner of fluid
		sequence = cell2struct( circshift( struct2cell(sequence),-1 ), stringSet(1:opt.nColumn) );

%       The simulation of columns within a SMB unit by the sequence, 
%       say, 'a', 'b', 'c', 'd' in four-column cases
        for j = 1:opt.nInterval

            for k = 1:opt.nColumn

                column = SMB.massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
                [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, column.initialState, varargin{:});

                currentData{eval(['sequence' '.' string(k)])}.outlet     = outletProfile;
                currentData{eval(['sequence' '.' string(k)])}.lastState  = lastState;

                tempData{eval(['sequence' '.' string(k)])}.concentration{j} = outletProfile.concentration;
            end

            if opt.nZone == 4
                dyncData{1, i} = currentData{ eval(['sequence' '.' char(stringSet(sum(opt.structID(1:3))))]) }.outlet.concentration;
                dyncData{2, i} = currentData{ eval( ['sequence' '.' char(stringSet(opt.structID(1)))]) }.outlet.concentration;
            elseif opt.nZone == 5
                dyncData{1, i} = currentData{ eval(['sequence' '.' char(stringSet(sum(opt.structID(1:4))))]) }.outlet.concentration;
                dyncData{2, i} = currentData{ eval(['sequence' '.' char(stringSet(sum(opt.structID(1:2))))]) }.outlet.concentration;
                dyncData{3, i} = currentData{ eval(['sequence' '.' char(stringSet(opt.structID(1)))])}.outlet.concentration;
            end

        end


%       Store the data of one round (opt.nColumn switches), into plotData
        for ii = 1:opt.nColumn
            transitionData{ii}.outlet.concentration = cat( 1, tempData{ii}.concentration{:} );
        end

        index = mod(i, opt.nColumn);
        if index == 0
            plotData(:,opt.nColumn) = transitionData';
        else
            plotData(:,index) = transitionData';
        end

%       Convergence criterion was adopted in each nColumn iteration
%           ||( C(z, t) - C(z, t + nColumn * t_s) ) / C(z, t)|| < tol, for a specific column
        if fix(i/opt.nColumn) == i/(opt.nColumn)

            diffNorm = 0; stateNorm = 0;

            for k = 1:opt.nComponents
                diffNorm = diffNorm + norm( convergPrevious(:,k) - transitionData{convergIndx}.outlet.concentration(:,k) );
                stateNorm = stateNorm + norm( transitionData{convergIndx}.outlet.concentration(:,k) );
            end

            relativeDelta = diffNorm / stateNorm;

            if opt.enableDebug
                fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
            end

%           Plot the outlet profile of each column in one round
            SMB.plotFigures(opt, plotData);
%           Plot the dynamic trajectory
            SMB.plotDynamic(opt, dyncData(:,1:i), i);

            if relativeDelta <= opt.tolIter
                break
            else
                convergPrevious = transitionData{convergIndx}.outlet.concentration;
            end

        end
    end
%-----------------------------------------------------------------------------------------


%   Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(plotData, opt);

%   Construct your own Objective Function and calculate the value
    objective = SMB.objectiveFunction(Results, opt);

    tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end

%   store the final data into DATA.mat file in the mode of forward simulation
    if opt.enableDebug
        save(sprintf('DATA_%2d.mat',fix(rand*100)),'Results');
        fprintf('The results have been stored in the DATA.mat \n');
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
% 
%      Copyright © 2008-2016: Eric von Lieres, Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================