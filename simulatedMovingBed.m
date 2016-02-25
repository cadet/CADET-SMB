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
    currentData  = cell(1, opt.nColumn);
	columnNumber = cell(1, opt.nColumn);

    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents);
    end

%   Number the columns for the sake of plotting 
    for k = 1:opt.nColumn
		columnNumber{k} = k;
    end

%	Specify the column for convergence checking	
    if opt.nZone == 4
        convergIndx = sum(opt.structID(1:2)) + 1;
    elseif opt.nZone == 5
        convergIndx = sum(opt.structID(1:3)) + 1;
    end

%   Preallocation
%   convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;

% 	The data for dynamic trajectory plotting
    if opt.nZone == 4
		dyncData = cell(2, opt.nMaxIter);
	elseif opt.nZone == 5
		dyncData = cell(3, opt.nMaxIter);
    end

% 	Construct the string in order to tell simulator the calculation sequence
	stringSet = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm'...
                 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'...
                 'aa' 'bb' 'cc' 'dd' 'ee' 'ff' 'gg' 'hh' 'ii' 'jj' 'kk' 'll' 'mm'...
                 'nn' 'oo' 'pp' 'qq' 'rr' 'ss' 'tt' 'uu' 'vv' 'ww' 'xx' 'yy' 'zz'};

	if opt.nColumn > length(stringSet)
		error('The simulation of %3g-column case in %3g-zone is not finished so far', opt.nColumn, opt.nZone);
	else
		string = stringSet(1:opt.nColumn);
	end

	sequence = cell2struct( columnNumber, string, 2 );

	string = char( [fliplr( string(1:convergIndx) ) fliplr( string(convergIndx+1:end) )] );


%----------------------------------------------------------------------------------------- 
%   Main loop
    for i = 1:opt.nMaxIter

        if mod(i, opt.nColumn) ~= 0
            k = mod(i, opt.nColumn);
        else
            k = opt.nColumn;
        end

        if i == 1
            initialState = [];
        end

%       The simulation of columns within a SMB unit by the sequence, 
%       say, 'c', 'b', 'a', 'd' in four-column cases
        column = SMB.massConservation(currentData, interstVelocity, Feed, opt, sequence, string(k));
        [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, initialState, varargin{:});

        currentData{eval(['sequence' '.' string(k)])}.outlet = outletProfile;
        initialState  = lastState;

        if mod(i, opt.nColumn) == 0
            if opt.nZone == 4
                dyncData{1, i/opt.nColumn} = currentData{eval(['sequence' '.' char(stringSet(sum(opt.structID(1:3))))])}.outlet.concentration;
                dyncData{2, i/opt.nColumn} = currentData{eval(['sequence' '.' char(stringSet(opt.structID(1)))])}.outlet.concentration;
            elseif opt.nZone == 5
                dyncData{1, i/opt.nColumn} = currentData{eval(['sequence' '.' char(stringSet(sum(opt.structID(1:4))))])}.outlet.concentration;
                dyncData{2, i/opt.nColumn} = currentData{eval(['sequence' '.' char(stringSet(sum(opt.structID(1:2))))])}.outlet.concentration;
                dyncData{3, i/opt.nColumn} = currentData{eval(['sequence' '.' char(stringSet(opt.structID(1)))])}.outlet.concentration;
            end
        end


%       Convergence criterion was adopted in each nColumn iteration
%           ||( C(z, t) - C(z, t + nColumn * t_s) ) / C(z, t)|| < tol, for a specific column
        if fix(i/opt.nColumn) == i/(opt.nColumn)

			diffNorm = 0; stateNorm = 0;

			for k = 1:opt.nComponents
				diffNorm = diffNorm + norm( convergPrevious(:,k) - currentData{convergIndx}.outlet.concentration(:,k) );
				stateNorm = stateNorm + norm( currentData{convergIndx}.outlet.concentration(:,k) );
			end

            relativeDelta = diffNorm / stateNorm;

            if opt.enableDebug
                fprintf('---- Round: %3d    Switch: %4d    CSS_relError: %g \n', i/opt.nColumn, i, relativeDelta);
            end

%           Plot the outlet profile of each column in one round
            SMB.plotFigures(opt, currentData);
%           Plot the dynamic trajectory
            SMB.plotDynamic(opt, dyncData(:,1:i/opt.nColumn), i/opt.nColumn);

            if relativeDelta <= opt.tolIter
                break
            else
                convergPrevious = currentData{convergIndx}.outlet.concentration;
            end

        end
    end
%-----------------------------------------------------------------------------------------


%   Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(currentData, opt);

%   Construct your own Objective Function and calculate the value
    objective = SMB.objectiveFunction(Results, opt);

    tTotal = toc(tTotal);
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
    end

%   Store the final data into DATA.mat file in the mode of forward simulation
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