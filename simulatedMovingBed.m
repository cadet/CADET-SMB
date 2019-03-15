function objective = simulatedMovingBed(varargin)
% =============================================================================
% Main function, which is in charge of switching to achieve CSS
%
% Binary separations (conventional four-zone scheme) and ternary separations
% (cascade, five-zone or eight zone schemes) are available in this program.
% As for the specific column configurations, I refer your to our paper.
% =============================================================================


    global string stringSet dummyProfile startingPointIndex Feed2;

    tTotal = tic;

    % Generate alphabet strings to name columns
    stringSet = SMB.stringGeneration();

    % Read operating parameters of SMB unit
    [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin{:});
    % If SMA isotherm, salt is considered as a component
    if strcmp('StericMassActionBinding', opt.BindingModel), opt.nComponents = opt.nComponents + 1; end

    % Check interstitial velocities in SMB optimization
    flag = SMB.interstVelocityCheck(interstVelocity, opt);
    % if flag is true (anyone is negative), assign a big objective value to get rid of
    if flag == 1, objective = 1e5; return; end

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
    string = char(stringSet(1:opt.nColumn));

    % Simulations follow the string sequence (starting from desorbent node)
    % To change starting simulation position, shift num to corresponding value
    string = circshift(string, 0);
    % Be aware, in eight-zone case, simulations cannot begain after feed2 node
    startingPointIndex = SMB.nodeIndexing(opt, string(1));

    % currentData matrix is the profile matrix of columns within one time interval
    currentData    = cell(1, opt.nColumn);
    % temp memory for transitionData
    tempData       = cell(1, opt.nColumn);
    % transtionData is the instant concentration profile of all columns in one iteration
    transitionData = cell(1, opt.nColumn);

    for k = 1:opt.nColumn
        currentData{k}.outlet.time = linspace(0, opt.switch, opt.timePoints);
        currentData{k}.outlet.concentration = zeros(length(Feed.time), opt.nComponents);
        currentData{k}.lastState = cell(1,2);

        if opt.enable_DPFR
            currentData{k}.lastState_DPFR = cell(1,2);
            currentData{k}.lastState_DPFR{1} = zeros(opt.nComponents, opt.DPFR_nCells); % DPFR before
            currentData{k}.lastState_DPFR{2} = zeros(opt.nComponents, opt.DPFR_nCells); % DPFR after
        end

        % This is the temporary data for dynamic trajectory plotting
        tempData{k}.concentration = cell(1, opt.nInterval);

        % transitionData stores the composition profile of nInterval pieces
        transitionData{k}.outlet.time = linspace(0, opt.switch*opt.nInterval, opt.timePoints*opt.nInterval);
        transitionData{k}.outlet.concentration = zeros(opt.timePoints*opt.nInterval, opt.nComponents);
    end

    % Generate arabic numbers to identify columns
    columnNumber = cell(1, opt.nColumn);
    for k = 1:opt.nColumn
        if k == 1, columnNumber{1} = opt.nColumn; else columnNumber{k} = k-1; end
    end
    % Combine alphabet string with arabic numbers for switching sake
    sequence = cell2struct( columnNumber, stringSet(1:opt.nColumn), 2 );

    % Specify the column for the convergence checking. The column after the Feed is usually adopted
    if opt.nZone == 4
        convergIndx = sum(opt.structID(1:2));
    elseif opt.nZone == 5
        convergIndx = sum(opt.structID(1:3));
    elseif opt.nZone == 8
        convergIndx = sum(opt.structID(1:3));
    else
        error('Please choose the correct zone configuration \n');
    end

    % convergPrevious is used for stopping criterion
    convergPrevious = transitionData{convergIndx}.outlet.concentration;
    % The profile of last column in terms of sequence is stored as dummyProfile
    dummyProfile    = currentData{sequence.(string(end))}.outlet;

    % plotData (columnNumber x switches), monitoring instant profile of each column in one iteration
    plotData = cell(opt.nColumn,opt.nColumn);

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

        % Switch implementation by means of attaching different column to corresponding position
        sequence = cell2struct( circshift( struct2cell(sequence),-1 ), stringSet(1:opt.nColumn) );

        simMex = cell(1, opt.nColumn);

        for j = 1:opt.nInterval

            for k = string' % do nColumn simulations in terms of string sequence

                % The node balance: transmission of concentration, column state, velocity and so on
                column = SMB.massConservation(currentData, interstVelocity, Feed, Desorbent, opt, sequence, k, j);

                if opt.enable_CSTR

                    % The CSTR before the current column
                    column.inlet = SMB.CSTR(column.inlet, column, opt);

                    if j == 1
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, column.initialState, sequence.(k), j, varargin{:});
                    elseif j == opt.nInterval
                        [outletProfile, simMex, lastState] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    else
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    end

                    % The CSTR after the current column
                    outletProfile.outlet = SMB.CSTR(outletProfile.outlet, column, opt);

                elseif opt.enable_DPFR

                    % The DPFR before the current column
                    [column.inlet, lastState_DPFR_pre] = SMB.DPFR(column.inlet, column.initialState_DPFR{1}, opt);

                    if j == 1
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, column.initialState, sequence.(k), j, varargin{:});
                    elseif j == opt.nInterval
                        [outletProfile, simMex, lastState] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    else
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    end

                    % The DPFR after the current column
                    [outletProfile.outlet, lastState_DPFR_pos] = SMB.DPFR(outletProfile.outlet, column.initialState_DPFR{2}, opt);

                    currentData{sequence.(k)}.lastState_DPFR = [{lastState_DPFR_pre}, {lastState_DPFR_pos}];

                else

                    % The simulation of a single column with the CADET solver
                    if j == 1
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, column.initialState, sequence.(k), j, varargin{:});
                    elseif j == opt.nInterval
                        [outletProfile, simMex, lastState] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    else
                        [outletProfile, simMex] = SMB.secColumn(column.inlet, column.params, simMex, [], sequence.(k), j, varargin{:});
                    end


                end

                % when j equals nInterval, the column state is complete
                if j == opt.nInterval
                    currentData{sequence.(k)}.lastState = lastState;
                    currentData{sequence.(k)}.colState  = outletProfile.column;
                end
                currentData{sequence.(k)}.outlet        = outletProfile.outlet;

                % Store the interval-wise concentration profile of all columns
                tempData{sequence.(k)}.concentration{j} = outletProfile.outlet.concentration;

                % If eight-zone case, store the intermediate outlet as the Feed2 inlet
                if opt.nZone == 8
                    if ( strcmp('raffinate', opt.intermediate_feed) && strcmp(k, stringSet(sum(opt.structID(1:3)))) ) ...
                            || ( strcmp('extract',opt.intermediate_feed) && strcmp(k, stringSet(opt.structID(1))) )
                        Feed2 = outletProfile.outlet;
                    end
                    % desalting step in SMA isotherm during connection
                    if strcmp('StericMassActionBinding', opt.BindingModel)
                        Feed2.concentration = Feed2.concentration .* interstVelocity.raffinate1 ./ interstVelocity.feed2;
                    end
                end

            end % for k = string'


            % The collection of the dyncData for the trajectory plotting
            if opt.enableDebug
                if opt.nZone == 4
                    dyncData{1, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3))))) }.outlet.concentration;
                    dyncData{2, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(opt.structID(1)))) }.outlet.concentration;
                elseif opt.nZone == 5
                    dyncData{1, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:4))))) }.outlet.concentration;
                    if strcmp('extract', opt.intermediate) % two extract ports
                        dyncData{2, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:2))))) }.outlet.concentration;
                    elseif strcmp('raffinate', opt.intermediate) % two raffinate ports
                        dyncData{2, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3))))) }.outlet.concentration;
                    end
                    dyncData{3, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(opt.structID(1)))) }.outlet.concentration;
                elseif opt.nZone == 8
                    dyncData{1, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:7))))) }.outlet.concentration;
                    dyncData{2, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:5))))) }.outlet.concentration;
                    if strcmp('extract', opt.intermediate_feed) % two extract ports
                        dyncData{3, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(sum(opt.structID(1:3))))) }.outlet.concentration;
                    elseif strcmp('raffinate', opt.intermediate_feed) % two raffinate ports
                        dyncData{3, j+(i-1)*opt.nInterval} = currentData{sequence.(char(stringSet(opt.structID(1)))) }.outlet.concentration;
                    end
                end
            end

            % Interactive plottings when debug mode is active
            SMB.UIplot('interval', opt, i, j);

            % Plot the dynamic trajectory
            if opt.enableDebug, SMB.plotDynamic(opt, dyncData(:,1:j+(i-1)*opt.nInterval), j+(i-1)*opt.nInterval); end

        end % for j = 1:opt.nInterval


        % Transfer the dummyProfile, which is the profile of last interval of last column
        dummyProfile = outletProfile.outlet;

        % Plot the outlet profile of each column in nColumn switches
        if opt.enableDebug, SMB.plotFigures(opt, currentData); end

        % Cat nInterval pieces of tempData into transitionData,
        % which is the internal column state after nColumn switches
        for ii = 1:opt.nColumn
            transitionData{ii}.outlet.concentration = cat( 1, tempData{ii}.concentration{:} );
        end

        % Store the instant outlet profiles of all columns in one iteration
        index = mod(i, opt.nColumn);
        if index == 0
            plotData(:,opt.nColumn) = transitionData';
        else
            plotData(:,index) = transitionData';
        end

        % Convergence criterion adopted in each nColumn switches
        %   ||( C(z, t) - C(z, t + nColumn * t_s) ) / C(z, t)|| < tol, for a specific column
        if mod(i, opt.nColumn) == 0

            diffNorm = 0; stateNorm = 0;

            for k = 1:opt.nComponents
                diffNorm = diffNorm + norm( convergPrevious(:,k) - transitionData{convergIndx}.outlet.concentration(:,k) );
                stateNorm = stateNorm + norm( transitionData{convergIndx}.outlet.concentration(:,k) );
            end

            relativeDelta = diffNorm / stateNorm;

            % Interactive plotting when debug mode is active
            SMB.UIplot('cycle', opt, i, j, relativeDelta);

            if relativeDelta <= opt.tolIter
                break
            else
                convergPrevious = transitionData{convergIndx}.outlet.concentration;
            end

        end

    end % main loop


%   Post-process
%----------------------------------------------------------------------------------------
    % Compute the performance index, such Purity and Productivity
    Results = SMB.Purity_Productivity(plotData, varargin{:});

    % Construct your own Objective Function and calculate the value
    objective = SMB.objectiveFunction(Results, opt);

    tTotal = toc(tTotal);
    % Store the final data into DATA.mat file when debug is active
    if opt.enableDebug
        fprintf('The time elapsed for reaching the Cyclic Steady State: %g sec \n', tTotal);
        SMB.concDataConvertToASCII(currentData, opt);
        SMB.trajDataConvertToASCII(dyncData, opt);
        c = clock;
        save(sprintf('Performance_%d%d.mat',c(3), c(4)), 'Results');
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
