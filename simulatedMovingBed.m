function objective = simulatedMovingBed(varargin)
% =============================================================================
% Main function, which is in charge of switching to achieve CSS
%
% Binary separations (conventional four-zone scheme) and ternary separations
% (cascade, five-zone or eight zone schemes) are available in this program.
% As for the specific column configurations, I refer your to our paper.
% =============================================================================


    global string stringSet structID Feed2;

    tTotal = tic;

    % Generate alphabet strings to name columns
    stringSet = SMB.stringGeneration();

    % Read operating parameters of SMB unit
    [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin{:});

    % Check interstitial velocities in SMB optimization
    flag = SMB.interstVelocityCheck(interstVelocity, opt);
    % if flag is true (anyone is negative), assign a big objective value to get rid of
    if flag == 1, objective = 1e5; return; end

    if ~isfield(opt, 'structID')
        opt.structID = structID;
    end

    if opt.enable_CSTR && opt.enable_DPFR
        error('It is not allowed have both the CSTR and DPFR in the simulation \n');
    end

    if opt.nColumn > length(stringSet)
        error('The simulation of %3g-column case in %3g-zone is not finished so far \n', opt.nColumn, opt.nZone);
    end


%   Preallocation
%----------------------------------------------------------------------------------------
    % Generate an initial alphabet string set to identify positions of SMB unit.
    % The position after desorbent node is, by default, marked as "a" (first one)
    string = stringSet(1:opt.nColumn);

    % pre-allocate memory to currentData matrix, which is the profile matrix of columns
    currentData = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents);
    end
    initialState = cell(1,2);

    if opt.enable_DPFR
        initialState_DPFR = cell(1,2);
        initialState_DPFR{1} = zeros(opt.nComponents, opt.DPFR_nCells); % DPFR before
        initialState_DPFR{2} = zeros(opt.nComponents, opt.DPFR_nCells); % DPFR after
    end

    % Generate arabic numbers to identify pseudo columns
    columnNumber = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        columnNumber{k} = k;
    end
    % Combine alphabet string with arabic numbers for switching sake
	sequence = cell2struct( columnNumber, string, 2 );

    % Specify the column for the convergence checking. The column after the Feed is usually adopted
    if opt.nZone == 4
        convergIndx = sum(opt.structID(1:2)) + 1;
    elseif opt.nZone == 5
        convergIndx = sum(opt.structID(1:3)) + 1;
    elseif opt.nZone == 8
        convergIndx = sum(opt.structID(1:3)) + 1;
    else
        error('Please choose the correct zone configuration \n');
    end

    % convergPrevious is used for stopping criterion
    convergPrevious = currentData{convergIndx}.outlet.concentration;
    % The column after Feed node is selected as the analogous column in OCA
    string = char( [ fliplr(string(1:convergIndx)) fliplr(string(convergIndx+1:end)) ] );

    % dyncData is used for generating trajectories of withdrawn ports
    if opt.enableDebug
        if opt.nZone == 4
            dyncData = cell(2, opt.nMaxIter);
        elseif opt.nZone == 5
            dyncData = cell(3, opt.nMaxIter);
        elseif opt.nZone == 8
            dyncData = cell(3, opt.nMaxIter);
        end
    end


%   Simulations
%----------------------------------------------------------------------------------------
    % Interactive plotting when debug mode is enabled
    SMB.UIplot('head', opt);

    % Main loop
    for i = 1:opt.nMaxIter

        % Switch implementation by means of attaching the analog column to corresponding position
        if mod(i, opt.nColumn) ~= 0
            k = string(mod(i, opt.nColumn));
        else
            k = string(opt.nColumn);
        end

        % The node balance: transmission of concentration, column state, velocity and so on
        column = SMB.massConservation(currentData, interstVelocity, Feed, Desorbent, opt, sequence, k);

        if opt.enable_CSTR

            % The CSTR before the current column
            column.inlet = SMB.CSTR(column.inlet, column, opt);

            [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, initialState, varargin{:});

            % The CSTR after the current column
            outletProfile.outlet = SMB.CSTR(outletProfile.outlet, column, opt);

        elseif opt.enable_DPFR

            % The DPFR before the current column
            [column.inlet, lastState_DPFR_pre] = SMB.DPFR(column.inlet, initialState_DPFR{1}, opt);

            [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, initialState, varargin{:});

            % The DPFR after the current column
            [outletProfile.outlet, lastState_DPFR_pos] = SMB.DPFR(outletProfile.outlet, initialState_DPFR{2}, opt);

            initialState_DPFR = [{lastState_DPFR_pre}, {lastState_DPFR_pos}];

        else

            % The simulation of a single column with the CADET solver
            [outletProfile, lastState] = SMB.secColumn(column.inlet, column.params, initialState, varargin{:});

        end

        currentData{sequence.(k)}.outlet   = outletProfile.outlet;
        currentData{sequence.(k)}.colState = outletProfile.column;
        initialState = lastState;

        % If eight-zone case, store the intermediate outlet as the Feed2 inlet
        if opt.nZone == 8
            if ( strcmp('raffinate', opt.intermediate_feed) && strcmp(k, stringSet(sum(opt.structID(1:3)))) ) ...
                    || ( strcmp('extract', opt.intermediate_feed) && strcmp(k, stringSet(opt.structID(1))) )
                Feed2 = outletProfile.outlet;
            end
        end

        % The collection of the dyncData for the trajectory plotting
        if opt.enableDebug && mod(i, opt.nColumn) == 0
            if opt.nZone == 4
                dyncData{1, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3)))))}.outlet.concentration;
                dyncData{2, i/opt.nColumn} = currentData{sequence.(char(stringSet(opt.structID(1))))}.outlet.concentration;
            elseif opt.nZone == 5
                dyncData{1, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:4)))))}.outlet.concentration;
                if strcmp('extract', opt.intermediate) % two extract ports
                    dyncData{2, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:2)))))}.outlet.concentration;
                elseif strcmp('raffinate', opt.intermediate) % two raffinate ports
                    dyncData{2, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3)))))}.outlet.concentration;
                end
                dyncData{3, i/opt.nColumn} = currentData{sequence.(char(stringSet(opt.structID(1))))}.outlet.concentration;
            elseif opt.nZone == 8
                dyncData{1, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:7)))))}.outlet.concentration;
                dyncData{2, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:5)))))}.outlet.concentration;
                if strcmp('extract', opt.intermediate_feed) % two extract ports
                    dyncData{3, i/opt.nColumn} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3)))))}.outlet.concentration;
                elseif strcmp('raffinate', opt.intermediate_feed) % two raffinate ports
                    dyncData{3, i/opt.nColumn} = currentData{sequence.(char(stringSet(opt.structID(1))))}.outlet.concentration;
                end
            end
        end

        % Convergence criterion adopted in each nColumn switches
        %   ||( C(z, t) - C(z, t + nColumn * t_s) ) / C(z, t)|| < tol, for a specific column
        if mod(i, opt.nColumn) == 0

            diffNorm = 0; stateNorm = 0;

            for k = 1:opt.nComponents
                diffNorm = diffNorm + norm( convergPrevious(:,k) - currentData{convergIndx}.outlet.concentration(:,k) );
                stateNorm = stateNorm + norm( currentData{convergIndx}.outlet.concentration(:,k) );
            end

            relativeDelta = diffNorm / stateNorm;

            % Interactive plotting when debug mode is active
            SMB.UIplot('cycle', opt, i, relativeDelta);

            % Plot the outlet profile of each column in each switching period and dynamic trajectories at withdrawn ports
            if opt.enableDebug, SMB.plotFigures(opt, currentData); SMB.plotDynamic(opt, dyncData(:,1:i/opt.nColumn), i/opt.nColumn); end

            if relativeDelta <= opt.tolIter
                break
            else
                convergPrevious = currentData{convergIndx}.outlet.concentration;
            end

        end

    end % main loop


%   Post-process
%----------------------------------------------------------------------------------------
    % Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(currentData, opt);

    % Construct your own Objective Function and calculate the value
    try
        objective = SMB.objectiveFunction(Results, opt);
    catch e
        fprintf('%s\n', e.message);
        objective = 68106800;
    end

    tTotal = toc(tTotal);
    % Store the final data into DATA.mat file when debug is active
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
        SMB.concDataConvertToASCII(currentData, opt);
        SMB.trajDataConvertToASCII(dyncData, opt);
        c = clock;
        save(sprintf('Performance_%d%d.mat', c(3), c(4)), 'Results');
        fprintf('The results about concentration profiles and the trajectories have been stored \n');
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
%
%      Copyright © 2008-2017: Eric von Lieres, Qiaole He
%
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
