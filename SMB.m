
classdef SMB < handle
% ========================================================================================
% This is the class of the functions for simulated moving bed.
% ========================================================================================


    methods (Static = true, Access = 'public')


        function [outletProfile, simMex, lastState] = secColumn(inletProfile, params, simMex, lastState, num, interval, ParSwarm)
% ----------------------------------------------------------------------------------------
% Simulation of the single column
%
% Parameters:
%       - inletProfile. Inlet time and corresponding concentration profile
%       - params. Get parameters for simulation
%       - simMex. The cell that contains the MexSimulator for columns to resume
%       - lastState. The last STATE from previous simulation of next simulation
%       - num. The value is an alphabet which is used to indicate column
%       - interval. The indicator of the piece of operators
%       - ParSwarm. In the optimization situations, this is the vector which contains
%           the optimized decision variables
%
% Returns:
%       - outletProfile. outlet time and corresponding concentration profile
%       - simMex. The cell that contains the MexSimulator for columns
%       - lastState. Record the last STATE which is used as the boundary condition
% ----------------------------------------------------------------------------------------


            if nargin < 7
                ParSwarm = [];
                if nargin < 4
                    lastState = cell(1,2);
                    if nargin < 3
                        error('SMB.secColumn: There are no enough inputs for carrying out Simulator in CADET \n');
                    end
                end
            end

            if isempty(params.initialBulk) && isempty(params.initialSolid) && isempty(lastState{1})
                error('SMB.secColumn: There are no Initial / Boundary Conditions for the Simulator \n');
            end

            % Read operating parameters of the SMB unit
            [opt, ~, ~] = getParameters(ParSwarm);

            % Prepare the inlet profile to the CADET model
            inlet = PiecewiseCubicPolyProfile.fromUniformData(...
                inletProfile.time, inletProfile.concentration);

            % Nonnegative limitation
            inlet.constant(inlet.constant < 0) = 0;
            inlet.linear(inlet.linear < 0) = 0;
            inlet.quadratic(inlet.quadratic < 0) = 0;
            inlet.cubic(inlet.cubic < 0) = 0;

            % when simMex{num} is not empty, the sim.resume is used to continue simulations in operator-splitting
            if isempty(simMex{num})

                % General rate model
                mGrm = SingleGRM();

                % Discretization
                mGrm.nComponents    = opt.nComponents;
                % if SMA isotherm, salt component is regarded as a component
                if strcmp('StericMassActionBinding', opt.BindingModel)
                    mGrm.nComponents = opt.nComponents + 1;
                    opt.nComponents = opt.nComponents + 1;
                end
                mGrm.nCellsColumn   = opt.nCellsColumn;
                mGrm.nCellsParticle = opt.nCellsParticle;
                mGrm.nBoundStates   = ones(mGrm.nComponents, 1);

                % Initial conditions
                mGrm.initialBulk  = params.initialBulk;
                mGrm.initialSolid = params.initialSolid;
                if nargin >= 3 && ~isempty(lastState{1})
                    mGrm.initStateY    = lastState{1};
                    mGrm.initStateYdot = lastState{2};
                end

                % Transport
                mGrm.dispersionColumn          = params.dispersionColumn;
                mGrm.interstitialVelocity      = params.interstitialVelocity;
                mGrm.filmDiffusion             = opt.filmDiffusion;
                mGrm.diffusionParticle         = opt.diffusionParticle;
                mGrm.diffusionParticleSurface  = opt.diffusionParticleSurface;

                % Geometry
                mGrm.columnLength        = opt.columnLength;
                mGrm.particleRadius      = opt.particleRadius;
                mGrm.porosityColumn      = opt.porosityColumn;
                mGrm.porosityParticle    = opt.porosityParticle;

                % Adsorption
                if strcmp(opt.BindingModel, 'LinearBinding')

                    mLinear = LinearBinding();
                    mLinear.kineticBinding = false;

                    % Adsorption parameters
                    mLinear.kA   = opt.KA;
                    mLinear.kD   = opt.KD;

                    mGrm.bindingModel = mLinear;

                elseif strcmp(opt.BindingModel, 'MultiComponentLangmuirBinding')

                    mLangmuir = LangmuirBinding();
                    mLangmuir.kineticBinding = false;

                    mLangmuir.kA   = opt.KA;
                    mLangmuir.kD   = opt.KD;
                    mLangmuir.qMax = opt.QMAX;

                    mGrm.bindingModel = mLangmuir;

                elseif strcmp(opt.BindingModel, 'MultiComponentBiLangmuirBinding')

                    mBiLangmuir = BiLangmuirBinding();
                    mBiLangmuir.kineticBinding = false;

                    mBiLangmuir.kA1   = opt.KA(1);
                    mBiLangmuir.kD1   = opt.KD(1);
                    mBiLangmuir.qMax1 = opt.QMAX(1);
                    mBiLangmuir.kA2   = opt.KA(2);
                    mBiLangmuir.kD2   = opt.KD(2);
                    mBiLangmuir.qMax2 = opt.QMAX(2);

                    mGrm.bindingModel = mBiLangmuir;

                elseif strcmp(opt.BindingModel, 'StericMassActionBinding')

                    mSma = StericMassActionBinding();
                    mSma.kineticBinding = false;

                    mSma.lambda = opt.LAMBDA;
                    mSma.kA = opt.KA;
                    mSma.kD = opt.KD;
                    mSma.nu = opt.NU;
                    mSma.sigma = opt.SIGMA;

                    mGrm.bindingModel = mSma;

                else

                    error('%s: it is not available yet \n', opt.BindingModel);

                end

                mGrm.inlet = inlet;
                mGrm.returnSolutionBulk = true;

                % Create and configure simulator
                sim = Simulator.create();

                sim.solutionTimes       = inletProfile.time;
                sim.sectionTimes        = inlet.breaks;
                sim.sectionContinuity   = inlet.continuity;
                sim.model               = mGrm;

            else

                % Load the simMex to resume the simulations for operator-splitting
                sim = simMex{num};

                sim.model.inlet         = inlet;
                sim.solutionTimes       = inletProfile.time;
                sim.sectionTimes        = inlet.breaks;
                sim.sectionContinuity   = inlet.continuity;

            end

            sim.nThreads            = opt.nThreads;
            sim.returnLastState     = true;
            sim.returnLastStateSens = false;
            sim.absTol              = opt.ABSTOL;
            sim.initStepSize        = opt.INIT_STEP_SIZE;
            sim.maxSteps            = opt.MAX_STEPS;

            % Run the simulation
            try
                if isempty(simMex{num})
                    % Store the simMex for resume simulations
                    simMex{num} = sim;
                    result = sim.run();
                else
                    result = sim.resume();
                end
            catch e
                % Something went wrong
                error('CADET:simulationFailed', 'Check your settings and try again. \n%s',e.message);
            end

            % Extract the outletProfile
            outletProfile.outlet.time = result.solution.time;
            outletProfile.outlet.concentration = result.solution.outlet{1};
            outletProfile.column = SMB.extractColumnMatrix(result.solution.bulk{1}, opt);
            % When interval equals nInterval, the column state is complete
            if interval == opt.nInterval
                lastState{1} = result.solution.lastState;
                lastState{2} = result.solution.lastStateDot;
            end


        end % secColumn

        function column = massConservation(currentData, interstVelocity, Feed, Desorbent, opt, sequence, alphabet, j_interval)
% ----------------------------------------------------------------------------------------
% This is the function to calculate the concentration changes on each node.
%
% Parameters:
%       - currentData. Which includes each column's outlet concentration
%       (time-dependent), and the last state (which records every component's concentration
%        in bulk phase and stationary phase, and used as the initial state for the next simulation).
%       - interstVelocity. The interstitial velocity of each column
%       - Feed. The initialied injection
%       - opt. Options
%       - sequence. During switching, the structure used for storing the sequence of columns
%       - alphabet. It is a character. It tells this subroutine to calculate the specified column
%
% Returns:
%   Preparation for next column simulation
%       - column.inlet. The new inlet concentration of each column, which is
%       obtained from mass conservation on each node.
%       - column.lastState.
%       - column.params. Set the initial Mobile and Solid concentration to the
%       Simulator (if there is no lastState given), and also store the interstitial velocity.
% ----------------------------------------------------------------------------------------


            global stringSet dummyProfile startingPointIndex Feed2;

            if nargin < 8
                error('SMB.massConservation: There are no enough arguments \n');
            end

            % Time points
            column.inlet.time = linspace(0, opt.switch, opt.timePoints);

            % Interpret current alphabet to node position of the SMB unit
            index = SMB.nodeIndexing(opt, alphabet);

            % Specify previous column to obtain inlet concentration
            if ~strcmp(alphabet, 'a')
                pre_alphabet = char(alphabet - 1);
            else
                pre_alphabet = char(stringSet(opt.nColumn));
            end

            idx_i = sequence.(alphabet);     % the number of current column
            idx_j = sequence.(pre_alphabet); % the number of the column before

            % Get the interstitial velocity and boundary conditions
            params = SMB.getParams(sequence, interstVelocity, opt, index, alphabet, pre_alphabet);
            % Update the intersitial velocity, boundary conditions
            column.params = params{idx_i};
            column.initialState = currentData{idx_i}.lastState;
            % If DPFRs were implemented, transfer the boundary conditions between switches
            if opt.enable_DPFR
                column.initialState_DPFR = currentData{idx_i}.lastState_DPFR;
            end

            % Calculate concentration of the column due to its position in the SMB unit
            switch index

                case 'D' % node DESORBENT

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Desorbent.concentration .* interstVelocity.desorbent) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Desorbent.concentration .* interstVelocity.desorbent) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                case 'D1' % node DESORBENT1

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Desorbent{1}.concentration .* interstVelocity.desorbent1) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Desorbent{1}.concentration .* interstVelocity.desorbent1) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                case 'D2' % node DESORBENT2

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Desorbent{2}.concentration .* interstVelocity.desorbent2) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Desorbent{2}.concentration .* interstVelocity.desorbent2) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                case 'F' % node FEED in four-zone and five-zone

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                case 'F1' % node FEED1 in eight-zone

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed1) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed1) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                case 'F2' % node FEED2 in eight-zone

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                            params{idx_j}.interstitialVelocity + Feed2.concentration .* interstVelocity.feed2) ...
                            ./ params{idx_i}.interstitialVelocity;
                    else
                        column.inlet.concentration = (dummyProfile.concentration{j_interval} .* ...
                            params{idx_j}.interstitialVelocity + Feed2.concentration .* interstVelocity.feed2) ...
                            ./ params{idx_i}.interstitialVelocity;
                    end

                otherwise % node EXTRACT; RAFFINATE; MIDDLE

                    %   C_i^in = C_{i-1}^out
                    if ~strcmp(startingPointIndex, index)
                        column.inlet.concentration = currentData{idx_j}.outlet.concentration;
                    else
                        column.inlet.concentration = dummyProfile.concentration{j_interval};
                    end

            end

        end % massConservation

        function params = getParams(sequence, interstVelocity, opt, index, alphabet, pre_alphabet)
% ----------------------------------------------------------------------------------------
% After each swtiching, the value of velocities and initial conditions are changed
%
% Parameters:
% 		- sequence. The alphabet of each zone and their column number
% 		- interstVelocity. The interstitial velocity of each zone
% 		- opt. option of parameterss
% 		- index. The capital letter used to indicate which zone current column situated in
% 		- alphabet. The letter of current calculated column, i.
% 		- pre_alphabet. The letter of column before the current calculated one, i-1.
%
% Returns:
% 		- params. interstitial velocity, axial dispersion and boundary conditions
% ----------------------------------------------------------------------------------------


            if nargin < 6
                error('SMB.getParams: There are no enough arguments \n');
            end

            if length(opt.dispersionColumn) ~= opt.nZone
                error('SMB.getParams: The dimension of dispersionColumn in getParameters routine is not correct \n');
            end

            params = cell(1, opt.nColumn);
            % Set the initial conditions to the solver, but when lastState is used, this setup will be ignored
            params{sequence.(alphabet)} = struct('initialBulk', zeros(1,opt.nComponents), 'initialSolid',...
                zeros(1,opt.nComponents), 'interstitialVelocity', [], 'dispersionColumn', []);
            params{sequence.(pre_alphabet)} = struct('initialBulk', zeros(1,opt.nComponents), 'initialSolid',...
                zeros(1,opt.nComponents), 'interstitialVelocity', [], 'dispersionColumn', []);

            % For steric mass action isotherm, as salt concentration is considered
            if strcmp('StericMassActionBinding', opt.BindingModel)
                params{sequence.(alphabet)}.initialSolid = [opt.LAMBDA, zeros(1, opt.nComponents-1)];

                % step-wise salt concentration in ion-exchange SMB
                switch index
                    case {'F' 'M_F' 'R' 'M_R' 'F1' 'M_F1' 'R1' 'M_R1'}
                        params{sequence.(alphabet)}.initialBulk = [opt.concentrationSalt(1,1), zeros(1,opt.nComponents-1)];
                    case {'D' 'M_D' 'E' 'M_E' 'D1' 'M_D1' 'E1' 'M_E1'}
                        params{sequence.(alphabet)}.initialBulk = [opt.concentrationSalt(1,2), zeros(1,opt.nComponents-1)];
                    case {'F2' 'M_F2' 'R2' 'M_R2'}
                        params{sequence.(alphabet)}.initialBulk = [opt.concentrationSalt(2,1), zeros(1,opt.nComponents-1)];
                    case {'D2' 'M_D2' 'E2' 'M_E2'}
                        params{sequence.(alphabet)}.initialBulk = [opt.concentrationSalt(2,2), zeros(1,opt.nComponents-1)];
                end
            end

            if opt.nZone == 4

                switch index

                    case {'D' 'M_D'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(1);
                    case {'E' 'M_E'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(2);
                    case {'F' 'M_F'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract + interstVelocity.feed;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(3);
                    case {'R' 'M_R'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract + interstVelocity.feed;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(4);

                end

            elseif opt.nZone == 5

                if strcmp('extract', opt.intermediate)

                    switch index

                        case {'D' 'M_D'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(1);
                        case {'E1' 'M_E1'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(2);
                        case {'E2' 'M_E2'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.extract1 - interstVelocity.extract2;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(3);
                        case {'F' 'M_F'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.desorbent + interstVelocity.raffinate;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.extract1 - interstVelocity.extract2;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(4);
                        case {'R' 'M_R'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.desorbent + interstVelocity.raffinate;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(5);

                    end

                elseif strcmp('raffinate', opt.intermediate)

                    switch index

                        case {'D' 'M_D'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(1);
                        case {'E' 'M_E'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(2);
                        case {'F' 'M_F'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.extract + interstVelocity.feed;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(3);
                        case {'R1' 'M_R1'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.desorbent + interstVelocity.raffinate2;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.extract + interstVelocity.feed;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(4);
                        case {'R2' 'M_R2'}
                            params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                            params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                                - interstVelocity.desorbent + interstVelocity.raffinate2;
                            params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(5);

                    end

                end

            elseif opt.nZone == 8

                switch index

                    case {'D1' 'M_D1'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent1;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(1);
                    case {'E1' 'M_E1'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(2);
                    case {'F1' 'M_F1'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(3);
                    case {'R1' 'M_R1'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(4);
                    case {'D2' 'M_D2'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1 + interstVelocity.desorbent2;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(5);
                    case {'E2' 'M_E2'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.desorbent1 + interstVelocity.raffinate2 - interstVelocity.feed2;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1 + interstVelocity.desorbent2;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(6);
                    case {'F2' 'M_F2'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.desorbent1 + interstVelocity.raffinate2;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.desorbent1 + interstVelocity.raffinate2 - interstVelocity.feed2;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(7);
                    case {'R2' 'M_R2'}
                        params{sequence.(alphabet)}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent1;
                        params{sequence.(pre_alphabet)}.interstitialVelocity = interstVelocity.recycle ...
                            - interstVelocity.desorbent1 + interstVelocity.raffinate2;
                        params{sequence.(alphabet)}.dispersionColumn = opt.dispersionColumn(8);

                end

            end

        end % getParams

        function stringSet = stringGeneration()
% ----------------------------------------------------------------------------------------
% Generate strings to mark columns in SMB system
% ----------------------------------------------------------------------------------------


            stringSet = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm'...
                         'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'...
                         'aa' 'bb' 'cc' 'dd' 'ee' 'ff' 'gg' 'hh' 'ii' 'jj' 'kk' 'll' 'mm'...
                         'nn' 'oo' 'pp' 'qq' 'rr' 'ss' 'tt' 'uu' 'vv' 'ww' 'xx' 'yy' 'zz'};
%                         'a1' 'b1' 'c1' 'd1' 'e1' 'f1' 'g1' 'h1' 'i1' 'j1' 'k1' 'l1' 'm1'...
%                         'n1' 'o1' 'p1' 'q1' 'r1' 's1' 't1' 'u1' 'v1' 'w1' 'x1' 'y1' 'z1'...
%                         'a2' 'b2' 'c2' 'd2' 'e2' 'f2' 'g2' 'h2' 'i2' 'j2' 'k2' 'l2' 'm2'...
%                         'n2' 'o2' 'p2' 'q2' 'r2' 's2' 't2' 'u2' 'v2' 'w2' 'x2' 'y2' 'z2'...
%                         'a3' 'b3' 'c3' 'd3' 'e3' 'f3' 'g3' 'h3' 'i3' 'j3' 'k3' 'l3' 'm3'...
%                         'n3' 'o3' 'p3' 'q3' 'r3' 's3' 't3' 'u3' 'v3' 'w3' 'x3' 'y3' 'z3'...
%                         'a4' 'b4' 'c4' 'd4' 'e4' 'f4' 'g4' 'h4' 'i4' 'j4' 'k4' 'l4' 'm4'...
%                         'n4' 'o4' 'p4' 'q4' 'r4' 's4' 't4' 'u4' 'v4' 'w4' 'x4' 'y4' 'z4' ...
%                         'a5' 'b5' 'c5' 'd5' 'e5' 'f5' 'g5' 'h5' 'i5' 'j5' 'k5' 'l5' 'm5'...
%                         'n5' 'o5' 'p5' 'q5' 'r5' 's5' 't5' 'u5' 'v5' 'w5' 'x5' 'y5' 'z5' ...
%                         'a6' 'b6' 'c6' 'd6' 'e6' 'f6' 'g6' 'h6' 'i6' 'j6' 'k6' 'l6' 'm6'...
%                         'n6' 'o6' 'p6' 'q6' 'r6' 's6' 't6' 'u6' 'v6' 'w6' 'x6' 'y6' 'z6' ...
%                         'a7' 'b7' 'c7' 'd7' 'e7' 'f7' 'g7' 'h7' 'i7' 'j7' 'k7' 'l7' 'm7'...
%                         'n7' 'o7' 'p7' 'q7' 'r7' 's7' 't7' 'u7' 'v7' 'w7' 'x7' 'y7' 'z7' ...
%                         'a8' 'b8' 'c8' 'd8' 'e8' 'f8' 'g8' 'h8' 'i8' 'j8' 'k8' 'l8' 'm8'...
%                         'n8' 'o8' 'p8' 'q8' 'r8' 's8' 't8' 'u8' 'v8' 'w8' 'x8' 'y8' 'z8' ...
%                         'a9' 'b9' 'c9' 'd9' 'e9' 'f9' 'g9' 'h9' 'i9' 'j9' 'k9' 'l9' 'm9'...
%                         'n9' 'o9' 'p9' 'q9' 'r9' 's9' 't9' 'u9' 'v9' 'w9' 'x9' 'y9' 'z9' ...

        end % stringGeneration

        function index = nodeIndexing(opt, alphabet)
% ----------------------------------------------------------------------------------------
% This function interprets the alphabet of columns into the position of SMB unit
%
% Parameters:
% 		- opt. options involving the parameters for models
% 		- alphabet. The letter is used for switching
%
% Returns:
% 		- index. The indexing letter used for calculation of the mass conservation
% ----------------------------------------------------------------------------------------


            global stringSet;

            if nargin < 2
                error('SMB.nodeIndexing: There are no enough arguments \n');
            end

            string = stringSet(1:opt.nColumn);

            % Preallocation
            % stringBlock is used for storing the alphabet in each zone
            stringBlock = cell(1, opt.nZone);

            % Separate the string into nZone cells
            stringBlock{1} = string(1:opt.structID(1));
            for k = 2:opt.nZone
                stringBlock{k} = string( sum(opt.structID(1:k-1))+1 : sum(opt.structID(1:k)) );
            end

            % Assign each alphabet with the node indexing letter
            if opt.nZone == 4

                if any( strcmp(alphabet, stringBlock{1}) )
                    if strcmp(alphabet, stringBlock{1}(1)), index = 'D'; else index = 'M_D'; end
                elseif any( strcmp(alphabet, stringBlock{2}) )
                    if strcmp(alphabet, stringBlock{2}(1)), index = 'E'; else index = 'M_E'; end
                elseif any( strcmp(alphabet, stringBlock{3}) )
                    if strcmp(alphabet, stringBlock{3}(1)), index = 'F'; else index = 'M_F'; end
                elseif any( strcmp(alphabet, stringBlock{4}) )
                    if strcmp(alphabet, stringBlock{4}(1)), index = 'R'; else index = 'M_R'; end
                end

            elseif opt.nZone == 5

                if strcmp('extract', opt.intermediate)

                    % For two extract scheme
                    if any( strcmp(alphabet, stringBlock{1}) )
                        if strcmp(alphabet, stringBlock{1}(1)), index = 'D';  else index = 'M_D'; end
                    elseif any( strcmp(alphabet, stringBlock{2}) )
                        if strcmp(alphabet, stringBlock{2}(1)), index = 'E1'; else index = 'M_E1'; end
                    elseif any( strcmp(alphabet, stringBlock{3}) )
                        if strcmp(alphabet, stringBlock{3}(1)), index = 'E2'; else index = 'M_E2'; end
                    elseif any( strcmp(alphabet, stringBlock{4}) )
                        if strcmp(alphabet, stringBlock{4}(1)), index = 'F';  else index = 'M_F'; end
                    elseif any( strcmp(alphabet, stringBlock{5}) )
                        if strcmp(alphabet, stringBlock{5}(1)), index = 'R';  else index = 'M_R'; end
                    end

                elseif strcmp('raffinate', opt.intermediate)

                    % For two raffinate scheme
                    if any( strcmp(alphabet, stringBlock{1}) )
                        if strcmp(alphabet, stringBlock{1}(1)), index = 'D';  else index = 'M_D'; end
                    elseif any( strcmp(alphabet, stringBlock{2}) )
                        if strcmp(alphabet, stringBlock{2}(1)), index = 'E';  else index = 'M_E'; end
                    elseif any( strcmp(alphabet, stringBlock{3}) )
                        if strcmp(alphabet, stringBlock{3}(1)), index = 'F';  else index = 'M_F'; end
                    elseif any( strcmp(alphabet, stringBlock{4}) )
                        if strcmp(alphabet, stringBlock{4}(1)), index = 'R1'; else index = 'M_R1'; end
                    elseif any( strcmp(alphabet, stringBlock{5}) )
                        if strcmp(alphabet, stringBlock{5}(1)), index = 'R2'; else index = 'M_R2'; end
                    end

                end

            elseif opt.nZone == 8

                if any( strcmp(alphabet, stringBlock{1}) )
                    if strcmp(alphabet, stringBlock{1}(1)), index = 'D1'; else index = 'M_D1'; end
                elseif any( strcmp(alphabet, stringBlock{2}) )
                    if strcmp(alphabet, stringBlock{2}(1)), index = 'E1'; else index = 'M_E1'; end
                elseif any( strcmp(alphabet, stringBlock{3}) )
                    if strcmp(alphabet, stringBlock{3}(1)), index = 'F1'; else index = 'M_F1'; end
                elseif any( strcmp(alphabet, stringBlock{4}) )
                    if strcmp(alphabet, stringBlock{4}(1)), index = 'R1'; else index = 'M_R1'; end
                elseif any( strcmp(alphabet, stringBlock{5}) )
                    if strcmp(alphabet, stringBlock{5}(1)), index = 'D2'; else index = 'M_D2'; end
                elseif any( strcmp(alphabet, stringBlock{6}) )
                    if strcmp(alphabet, stringBlock{6}(1)), index = 'E2'; else index = 'M_E2'; end
                elseif any( strcmp(alphabet, stringBlock{7}) )
                    if strcmp(alphabet, stringBlock{7}(1)), index = 'F2'; else index = 'M_F2'; end
                elseif any( strcmp(alphabet, stringBlock{8}) )
                    if strcmp(alphabet, stringBlock{8}(1)), index = 'R2'; else index = 'M_R2'; end
                end

            end % opt.nZone

        end % nodeIndexing

        function obj = positionIndexing(opt)
% ----------------------------------------------------------------------------------------
% This is the function that tells the Purity_Productivity which column is used for
% calculation of purity and productivity. It is vital in arbitrary column
% configurations and it also depends on the plotData and the referred column.
%
% Parameters:
% 		- opt. options
%
% Returns:
% 		- obj. A struct data which contains the numbers for indicating the EXTRACT and
% 			RAFFINATE ports those are used for calculation of purity and productivity
% ----------------------------------------------------------------------------------------


            Number = circshift( (fliplr(1:opt.nColumn))', 1 );

            % numberBlock is used for storing the column number in each zone
            numberBlock = cell(1, opt.nZone);

            % Separate the number string into nZone cells
            numberBlock{1} = Number(1:opt.structID(1));
            for k = 2:opt.nZone
                numberBlock{k} = Number( sum(opt.structID(1:k-1))+1 : sum(opt.structID(1:k)) );
            end

            if opt.nZone == 4

                obj.position = [numberBlock{1}(end), numberBlock{3}(end)]; % [Desorbent, Feed] node

            elseif opt.nZone == 5

                if strcmp('extract', opt.intermediate)
                    obj.position = [numberBlock{1}(end), numberBlock{2}(end), numberBlock{4}(end)]; % [Desorbent, Extract1, Feed] node
                elseif strcmp('raffinate', opt.intermediate)
                    obj.position = [numberBlock{1}(end), numberBlock{3}(end), numberBlock{4}(end)]; % [Desorbent, Feed, Raffinate1] node
                end

            elseif opt.nZone == 8

                if strcmp('raffinate', opt.intermediate_feed)
                    obj.position = [numberBlock{1}(end), numberBlock{5}(end), numberBlock{7}(end)]; % [Desorbent1, Desorbent2, Feed2] node
                elseif strcmp('extract', opt.intermediate_feed)
                    obj.position = [numberBlock{3}(end), numberBlock{5}(end), numberBlock{7}(end)]; % [Feed1, Desorbent2, Feed2] node
                end

            end

        end % positionIndexing

        function flag = interstVelocityCheck(interstVelocity, opt)
% ----------------------------------------------------------------------------------------
% This function checks the existance of the negative velocity and returns false to flag
%
% Parameters:
% 		- opt. options
%       - interstVelocity. The interstitial velocity in the SMB unit
%
% Returns:
% 		- flag. flag = 1, there is negative velocity in the struct of interstVelocity
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.interstVelocity: There are no enough arguments \n');
            end

            if opt.nZone == 4

                velocity = [interstVelocity.recycle, interstVelocity.feed ...
                    interstVelocity.raffinate, interstVelocity.desorbent, interstVelocity.extract];

                flag = any(velocity <= 0);
                flag = flag || (interstVelocity.recycle - interstVelocity.desorbent) < 0;
                flag = flag || (interstVelocity.recycle - interstVelocity.extract) < 0;

            elseif opt.nZone == 5

                if strcmp('extract', opt.intermediate)

                    velocity = [interstVelocity.recycle, interstVelocity.feed ...
                        interstVelocity.raffinate, interstVelocity.desorbent ...
                        interstVelocity.extract1, interstVelocity.extract2];

                    flag = any(velocity <= 0);
                    flag = flag || (interstVelocity.recycle - interstVelocity.desorbent) < 0;
                    flag = flag || (interstVelocity.recycle - interstVelocity.extract1) < 0;
                    flag = flag || (interstVelocity.recycle - interstVelocity.extract1 - interstVelocity.extract2) < 0;

                elseif strcmp('raffinate', opt.intermediate)

                    velocity = [interstVelocity.recycle, interstVelocity.feed ...
                        interstVelocity.raffinate1, interstVelocity.raffinate2 ...
                        interstVelocity.desorbent, interstVelocity.extract];

                    flag = any(velocity <= 0);
                    flag = flag || (interstVelocity.recycle - interstVelocity.desorbent) < 0;
                    flag = flag || (interstVelocity.recycle - interstVelocity.extract) < 0;
                    flag = flag || (interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed - interstVelocity.raffinate1) < 0;

                end

            elseif opt.nZone == 8

                velocity = [interstVelocity.recycle, interstVelocity.desorbent1 ...
                    interstVelocity.extract1, interstVelocity.feed1 ...
                    interstVelocity.raffinate1, interstVelocity.desorbent2 ...
                    interstVelocity.extract2, interstVelocity.raffinate2];

                flag = any(velocity <= 0);
                flag = flag || (interstVelocity.recycle - interstVelocity.desorbent1) < 0;
                flag = flag || (interstVelocity.recycle - interstVelocity.extract1) < 0;
                flag = flag || (interstVelocity.recycle - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1) < 0;
                flag = flag || (interstVelocity.recycle - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1 + interstVelocity.desorbent2 - interstVelocity.extract2) < 0;

            end

        end % interstVelocityCheck

        function intervalAmountCheck(opt, interstVelocity)
%-----------------------------------------------------------------------------------------
% This is the function that checks the amount of the time sections
%
% Parameters:
% 		- opt. options
%       - interstVelocity. The interstitial velocity in the SMB unit
%-----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.intervalAmountCheck: There are no enough arguments \n');
            end

            if opt.nZone == 4

                flowRate_zone_I = interstVelocity.recycle;
                flowRate_zone_II = interstVelocity.recycle - interstVelocity.extract;
                flowRate_zone_III = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
                flowRate_zone_IV = interstVelocity.recycle - interstVelocity.desorbent;

                velocityAverage = (flowRate_zone_I + flowRate_zone_II + flowRate_zone_III + flowRate_zone_IV) / opt.nZone;

            elseif opt.nZone == 5

                flowRate_zone_I = interstVelocity.recycle;
                flowRate_zone_II = interstVelocity.recycle - interstVelocity.extract1;
                flowRate_zone_III = interstVelocity.recycle - interstVelocity.extract1 - interstVelocity.extract2;
                flowRate_zone_IV = interstVelocity.recycle - interstVelocity.desorbent + interstVelocity.raffinate;
                flowRate_zone_V = interstVelocity.recycle - interstVelocity.desorbent;

                velocityAverage = (flowRate_zone_I + flowRate_zone_II + flowRate_zone_III + flowRate_zone_IV + flowRate_zone_V) / opt.nZone;

            elseif opt.nZone == 8

                flowRate_zone_I = interstVelocity.recycle;
                flowRate_zone_II = interstVelocity.recycle - interstVelocity.extract1;
                flowRate_zone_III = interstVelocity.recycle - interstVelocity.extract1 + interstVelocity.feed1;
                flowRate_zone_IV = interstVelocity.recycle - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1;
                flowRate_zone_V = interstVelocity.recycle - interstVelocity.extract1 + interstVelocity.feed1 - interstVelocity.raffinate1 + interstVelocity.desorbent2;
                flowRate_zone_VI = interstVelocity.recycle - interstVelocity.desorbent1 + interstVelocity.raffinate2 - interstVelocity.feed2;
                flowRate_zone_VII = interstVelocity.recycle - interstVelocity.desorbent1 + interstVelocity.raffinate2;
                flowRate_zone_VIII = interstVelocity.recycle - interstVelocity.desorbent1;

                velocityAverage = (flowRate_zone_I + flowRate_zone_II + flowRate_zone_III + flowRate_zone_IV + flowRate_zone_V +...
                    flowRate_zone_VI + flowRate_zone_VII + flowRate_zone_VIII) / opt.nZone;

            end

            dummyTime = opt.columnLength / velocityAverage / 2;

            minimalInterval = ceil(opt.switch * opt.nInterval / dummyTime);

            if opt.nInterval <= minimalInterval
                error('getParameters: Please set the opt.nInterval in the getParameters larger than %4d \n', minimalInterval);
            end

%            if mod(opt.timePoints, opt.nInterval) ~= 0
%                error('getParameters: Please make the timePoints be divisible to the nInterval \n');
%            end

        end % intervalAmountCheck


        function [colState, lastState] = observerSimulation(sequence, interstVelocity, Feed, Desorbent, dummyProfile, initialState, alphabet, opt)
%----------------------------------------------------------------------------------------
% dummy column simulation
% In order to obtain the actual column state of the very first column, a dummy simulation
% is implemented with the outlet profile from previous column and initial column state.
%
% Parameters:
% 		- dummyProfile. In order to get the inlet profiel
% 		- initialState. The initially column state
% 		- alphabet. The index used for indicating which zone it is
% 		- Feed. The feed concentration profile
% 		- sequence. The column configuration
% 		- interstVelocity. It is used for calculating the interstitial velocity in each zone
% 		- opt. The parameters that might be useful in calculation
%
% Return:
% 		- lastState. The updated column state, also actual state of the very first column
%----------------------------------------------------------------------------------------


            global stringSet Feed2;

            if nargin < 8
                error('SMB.observerSimulation: There are no enough arguments \n');
            end

            index = SMB.nodeIndexing(opt, alphabet);

            if ~strcmp(alphabet, 'a')
                pre_alphabet = char(alphabet - 1);
            else
                pre_alphabet = char(stringSet(opt.nColumn));
            end

            % Get the interstitial velocity of each column and boundary conditions
            params = SMB.getParams(sequence, interstVelocity, opt, index, alphabet, pre_alphabet);
            idx_i  = sequence.(alphabet);     % the number of current column
            idx_j  = sequence.(pre_alphabet); % the number of the column before

            switch index

                case 'D' % node DESORBENT

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Desorbent.concentration .* interstVelocity.desorbent) ...
                        ./ params{idx_i}.interstitialVelocity;

                case 'D1' % node DESORBENT1

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Desorbent{1}.concentration .* interstVelocity.desorbent1) ...
                        ./ params{idx_i}.interstitialVelocity;

                case 'D2' % node DESORBENT2

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Desorbent{2}.concentration .* interstVelocity.desorbent2) ...
                        ./ params{idx_i}.interstitialVelocity;

                case 'F' % node FEED of four-zone and five-zone

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                        ./ params{idx_i}.interstitialVelocity;

                case 'F1' % node FEED1 of eight-zone

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed1) ...
                        ./ params{idx_i}.interstitialVelocity;

                case 'F2' % node FEED2

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    dummyProfile.concentration = (dummyProfile.concentration .* ...
                        params{idx_j}.interstitialVelocity + Feed2.concentration .* interstVelocity.feed2) ...
                        ./ params{idx_i}.interstitialVelocity;

                otherwise % node EXTRACT; RAFFINATE; MIDDLE

                    %   C_i^in = C_{i-1}^out
                    dummyProfile.concentration = dummyProfile.concentration;

            end

            % Simulation of the first column again to accelarate the convergence
            [outletProfile, ~, lastState] = SMB.secColumn(dummyProfile, params{idx_i}, cell(1,opt.nColumn), initialState, 1, opt.nInterval);
            colState = outletProfile.column;


        end % observerSimulation

        function objective = objectiveFunction(Results, opt)
% ----------------------------------------------------------------------------------------
% The objective function for the optimizers
% You can also define your own objective function here. The default function is:
%
% Max Productivity_extract + Productivity_raffinate
% s.t. Purity_extract   >= 99% for more retained component
%      Purity_raffinate >= 99% for less retained component
%      other implicit constraints, such as upbound on Desorbent consumption
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.objectiveFunction: There are no enough arguments \n');
            end

            if opt.nZone == 4
                % Construct the Penalty Function for the objective function
                penalty = abs( min(Results.Purity_1 - opt.Purity_limit(1), 0) ) * 100 * opt.Penalty_factor ...
                        + abs( min(Results.Purity_2 - opt.Purity_limit(2), 0) ) * 100 * opt.Penalty_factor;

                % (-) since in the optimizer, the defined program is of optimization of minimum
                objective = -(Results.Productivity_1 + Results.Productivity_2) + penalty;

            elseif opt.nZone == 5 || opt.nZone == 8
                % Construct the Penalty Function for the objective function
                penalty = abs( min(Results.Purity_1 - opt.Purity_limit(1), 0) ) * 100 * opt.Penalty_factor ...
                        + abs( min(Results.Purity_2 - opt.Purity_limit(2), 0) ) * 100 * opt.Penalty_factor ...
                        + abs( min(Results.Purity_3 - opt.Purity_limit(3), 0) ) * 100 * opt.Penalty_factor;

                % (-) since in the optimizer, the defined program is of optimization of minimum
                objective = -(Results.Productivity_1 + Results.Productivity_2 + Results.Productivity_3) + penalty;

            end

            if opt.enableDebug, fprintf('**** The objective value:  %g \n', objective); end

        end % objectiveFunction

        function Results = Purity_Productivity(plotData, varargin)
% ----------------------------------------------------------------------------------------
% Calculation of the performance index of SMB, such Purity and Productivity
%
% Parameters:
% 		- plotData. The data for plotting
%
% Returns:
% 		- Resutls. A struct data which contains the purity and productivity of extract
% 			and raffinate ports, respectively
% ----------------------------------------------------------------------------------------


            if nargin < 1
                error('SMB.Purity_Productivity: There are no enough arguments \n');
            end

            opt = getParameters(varargin{:});

            profile = SMB.squeezePlotData(plotData, opt);

            % The nominator of the formualar of productivity
            Nominator = pi * (opt.columnDiameter/2)^2 * opt.columnLength * (1-opt.porosityColumn);

            if opt.nZone == 4

                % 1 stands for extract while 2 stands for raffinate
                sum_1 = 0; sum_2 = 0;
                for k = 1:opt.nComponents
                    sum_1 = sum_1 + trapz(profile{1}.time, profile{1}.concentration(:,k));
                    sum_2 = sum_2 + trapz(profile{2}.time, profile{2}.concentration(:,k));
                end

                % Extract and raffiante ports
                Purity_1 = trapz(profile{1}.time, profile{1}.concentration(:, profile{1}.id)) / sum_1;
                Purity_2 = trapz(profile{2}.time, profile{2}.concentration(:, profile{2}.id)) / sum_2;

                % per switch time in the extract port, such (unit: mol/m^3/s) amount of target component was collected
                Productivity_1 = trapz(profile{1}.time, profile{1}.concentration(:, profile{1}.id)) ...
                    * opt.molMass(profile{1}.id) * profile{1}.flowrate / Nominator;
                Productivity_2 = trapz(profile{2}.time, profile{2}.concentration(:, profile{2}.id)) ...
                    * opt.molMass(profile{2}.id) * profile{1}.flowrate / Nominator;

                if opt.enableDebug
                    fprintf('Purity of extract: %g %% \n', Purity_1 * 100);
                    fprintf('Purity of raffinate: %g %% \n', Purity_2 * 100)
                    fprintf('Productivity of extract in each switching time: %g mol/m^3/s \n', Productivity_1);
                    fprintf('Productivity of raffinate in each switching time: %g mol/m^3/s \n', Productivity_2);
                end

                Results = struct('Purity_1', Purity_1, 'Purity_2', Purity_2, ...
                    'Productivity_1', Productivity_1, 'Productivity_2', Productivity_2);


            elseif opt.nZone == 5 || opt.nZone == 8

                sum_1 = 0; sum_2 = 0; sum_3 = 0;
                for k = 1:opt.nComponents
                    sum_1 = sum_1 + trapz(profile{1}.time, profile{1}.concentration(:,k));
                    sum_2 = sum_2 + trapz(profile{2}.time, profile{2}.concentration(:,k));
                    sum_3 = sum_3 + trapz(profile{3}.time, profile{3}.concentration(:,k));
                end

                % Extract and raffinate ports
                Purity_1 = trapz(profile{1}.time, profile{1}.concentration(:, profile{1}.id)) / sum_1;
                Purity_2 = trapz(profile{2}.time, profile{2}.concentration(:, profile{2}.id)) / sum_2;
                Purity_3 = trapz(profile{3}.time, profile{3}.concentration(:, profile{3}.id)) / sum_3;

                % per switch time in the tank of extract port, such (unit: mol/m^3/s) amount of target component was collected
                Productivity_1 = trapz(profile{1}.time, profile{1}.concentration(:, profile{1}.id))...
                    * opt.molMass(profile{1}.id) * profile{1}.flowrate / Nominator;
                Productivity_2 = trapz(profile{2}.time, profile{2}.concentration(:, profile{2}.id))...
                    * opt.molMass(profile{2}.id) * profile{2}.flowrate / Nominator;
                Productivity_3 = trapz(profile{3}.time, profile{3}.concentration(:, profile{3}.id))...
                    * opt.molMass(profile{3}.id) * profile{3}.flowrate / Nominator;

                if opt.enableDebug
                    fprintf('Purity of outlet 1: %g %% \n', Purity_1 * 100);
                    fprintf('Purity of outlet 2: %g %% \n', Purity_2 * 100);
                    fprintf('Purity of outlet 3: %g %% \n', Purity_3 * 100)
                    fprintf('Productivity of outlet 1 in each switching time: %g mol/m^3/s \n', Productivity_1);
                    fprintf('Productivity of outlet 2 in each switching time: %g mol/m^3/s \n', Productivity_2);
                    fprintf('Productivity of outlet 3 in each switching time: %g mol/m^3/s \n', Productivity_3);
                end

                Results = struct('Purity_1', Purity_1, 'Purity_2', Purity_2,...
                    'Purity_3', Purity_3, 'Productivity_1', Productivity_1,...
                    'Productivity_2', Productivity_2, 'Productivity_3', Productivity_3);

            end % opt.nZone

        end % Purity_Productivity

        function profile = squeezePlotData(plotData, opt)
% ----------------------------------------------------------------------------------------
% Squeeze the plotData in order to be used by purity & productivity subroutine
%
% Parameters:
% 		- plotData. The data for plotting
% 		- opt. options of parameters
%
% Returns:
% 		- Resutls. A struct data which contains the purity and productivity of extract
% 			and raffinate ports, respectively
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.squeezePlotData: There are no enough arguments \n');
            end

            obj = SMB.positionIndexing(opt);

            if opt.nZone == 4

                profile = cell(1,2);

                % 1 standards for extract port while 2 stands for raffinate port
                profile{1} = plotData{1, obj.position(1)}.outlet;
                profile{2} = plotData{1, obj.position(2)}.outlet;

                profile{1}.id = opt.compTargID(1);
                profile{2}.id = opt.compTargID(2);

                profile{1}.flowrate = opt.flowRate_extract;
                profile{2}.flowrate = opt.flowRate_raffinate;

                if strcmp('StericMassActionBinding', opt.BindingModel)
                    profile{1}.concentration = profile{1}.concentration(:, 2:end);
                    profile{2}.concentration = profile{2}.concentration(:, 2:end);
                end

            elseif opt.nZone == 5

                profile = cell(1,3);

                profile{1} = plotData{1, obj.position(1)}.outlet;
                profile{2} = plotData{1, obj.position(2)}.outlet;
                profile{3} = plotData{1, obj.position(3)}.outlet;

                profile{1}.id = opt.compTargID(1);
                profile{2}.id = opt.compTargID(2);
                profile{3}.id = opt.compTargID(3);

                if strcmp('extract', opt.intermediate)
                    profile{1}.flowrate = opt.flowRate_extract1;
                    profile{2}.flowrate = opt.flowRate_extract2;
                    profile{3}.flowrate = opt.flowRate_raffinate;
                elseif strcmp('raffinate', opt.intermediate)
                    profile{1}.flowrate = opt.flowRate_extract;
                    profile{2}.flowrate = opt.flowRate_raffinate1;
                    profile{3}.flowrate = opt.flowRate_raffinate2;
                end

                if strcmp('StericMassActionBinding', opt.BindingModel)
                    profile{1}.concentration = profile{1}.concentration(:, 2:end);
                    profile{2}.concentration = profile{2}.concentration(:, 2:end);
                    profile{3}.concentration = profile{3}.concentration(:, 2:end);
                end

            elseif opt.nZone == 8

                profile = cell(1,3);

                profile{1} = plotData{1, obj.position(1)}.outlet;
                profile{2} = plotData{1, obj.position(2)}.outlet;
                profile{3} = plotData{1, obj.position(3)}.outlet;

                profile{1}.id = opt.compTargID(1);
                profile{2}.id = opt.compTargID(2);
                profile{3}.id = opt.compTargID(3);

                if strcmp('raffinate', opt.intermediate_feed)
                    profile{1}.flowrate = opt.flowRate_extract1;
                    profile{2}.flowrate = opt.flowRate_extract2;
                    profile{3}.flowrate = opt.flowRate_raffinate2;
                elseif strcmp('extract', opt.intermediate_feed)
                    profile{1}.flowrate = opt.flowRate_raffinate1;
                    profile{2}.flowrate = opt.flowRate_extract2;
                    profile{3}.flowrate = opt.flowRate_raffinate2;
                end

                if strcmp('StericMassActionBinding', opt.BindingModel)
                    profile{1}.concentration = profile{1}.concentration(:, 2:end);
                    profile{2}.concentration = profile{2}.concentration(:, 2:end);
                    profile{3}.concentration = profile{3}.concentration(:, 2:end);
                end

            end % opt.nZone

        end % squeezePlotData

        function array = extractColumnMatrix(colState, opt)
% ----------------------------------------------------------------------------------------
% Squeece 3D column matrix to a 2D matrix, whose dimension is nCellsColumn X nComponent
% ----------------------------------------------------------------------------------------


            array = squeeze(colState(end,:,:));

        end % extractColumnMatrix

        function columnSpline = CSTR(Profile, column, opt)
% ----------------------------------------------------------------------------------------
% Simulation of the continuous stirred tank reactor (CSTR)
%
% Parameters:
%       - Profile. Inlet time concentration, the initial conditions
%       - column. The boundary conditions of the CSTR
%       - opt. Options for the software
%
% Returns:
%       - columnProfile. outlet time and corresponding concentration profile
%               - time. The time points observed
%               - concentration. The concentration of the outlet of CSTR due to the time
%               points
% ----------------------------------------------------------------------------------------


            if nargin < 3
                error('SMB.CSTR: There are no enough arguments \n');
            end

            for i = 1:opt.nComponents

                IC = Profile.concentration(1,i);

                [T, Y] = ode15s(@(t, y) ode_CSTR(t, y, Profile.time, Profile.concentration, column.params.interstitialVelocity, i), [0, opt.switch], IC);

                columnSpline.concentration(:,i) = interp1(T, Y, Profile.time);

            end

            columnSpline.time = Profile.time;


            function DyDt = ode_CSTR(t, y, time, concentration, interstVelocity, comp)
            % ----------------------------------------------------------------------------
            % The ODE function of continuous stirred tank reactor
            % DyDt = [ y_{out}(t) - y_{in}(t) ] / tau
            %
            % Parameters:
            %       - t. Time point
            %       - y. Concentration
            %       - time. The vector of the observed time points
            %       - concentration. The initial conditions for the IVP
            %       - interstVelocity. The velocity of CSTR
            %       - comp. The index of components
            % ----------------------------------------------------------------------------


                concentration_f = interp1(time, concentration(:,comp), t);

                tau = opt.CSTR_length / interstVelocity;

                DyDt = 1/tau * (concentration_f - y);

            end

        end % CSTR

        function [columnSpline, lastState]  = DPFR(Profile, initialState, opt)
% ----------------------------------------------------------------------------------------
% Simulation of the dispersive plug flow reactor (DPFR)
%
% Parameters:
%       - Profile. Inlet time concentration, the initial conditions
%       - initialState. The boundary conditions of the DPFR
%       - opt. Options for the software
%
% Returns:
%       - columnProfile. outlet time and corresponding concentration profile
%               - time. The time points observed
%               - concentration. The concentration of the outlet of CSTR due to the time
%               points
%       - lastState. Record the last state of the DPFRs as the initial state of next
%       calculation
% ----------------------------------------------------------------------------------------


            if nargin < 3
                error('SMB.DPFR: There are no enough arguments \n');
            end

            lastState = zeros(opt.nComponents, opt.DPFR_nCells);

			delta_h = opt.DPFR_length / opt.DPFR_nCells;

            for i = 1:opt.nComponents

                concentration = interp1(Profile.time, Profile.concentration(:,i), 'pchip', 'pp');

                [T, Y] = ode15s(@ode_DPFR, Profile.time, initialState(i,:));

                columnSpline.concentration(:,i) = interp1(T, Y(:,end), Profile.time);

                lastState(i, :) = Y(end, :);

            end

            columnSpline.time = Profile.time;


            function y = ode_DPFR(t, x)
            % ----------------------------------------------------------------------------
            % The PDE function of dispersive plug flow reactor
            % Dy/Dt = - u Dy/Dz + D_{ax} DDy/DDz
            %
            % Parameters:
            %       - t. Time point
            %       - x. Initial conditons
            % ----------------------------------------------------------------------------


                y = zeros(opt.DPFR_nCells, 1);

                % Convection
                y(1) = y(1) - opt.DPFR_velocity / delta_h * (x(1) - ppval(concentration, t));
                y(2:end) = y(2:end) - opt.DPFR_velocity / delta_h * (x(2:end) - x(1:end-1));

                % Dispersion
                y(1) = y(1) + opt.DPFR_dispersion / sqrt(delta_h) * (x(2) - x(1));
                y(2:end-1) = y(2:end-1) + opt.DPFR_dispersion / sqrt(delta_h) * (x(3:end) - 2 * x(2:end-1) + x(1:end-2));
                y(end) = y(end) + opt.DPFR_dispersion / sqrt(delta_h) * (x(end-1) - x(end));

            end

        end % DPFR

        function concDataConvertToASCII(plotData, opt)
% ---------------------------------------------------------------------------------------
% This fucntion converts the struct data in Matlab into ASCII format
%
% Parameters:
%       - plotData. The data that stores the concentration profile information of all the
%       columns. After conversion, the plotData will be divided into nColumn separate *.dat files
%       - opt. The option that might be used
% ---------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.concDataConvertToASCII: There are no enough arguments \n');
            end

            y = plotData{1,1}.colState;
            for k = 2:opt.nColumn
                y = [y; plotData{1,k}.colState];
            end

            save('chromatogram.dat', 'y' ,'-ascii', '-double', '-tabs');

        end % concDataConvertToASCII

        function trajDataConvertToASCII(dyncData, opt)
% ----------------------------------------------------------------------------------------
% This is a fucntion that converts the struct data in Matlab into ASCII format
%
% Parameters:
%       - dyncData. The data that stores the trajectory informtion
%       - opt. The option that might be used
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.trajDataConvertToASCII: There are no enough arguments \n');
            end

            % Find the index of non-empty cell
            for i = 1:length(dyncData)
                if isempty(dyncData{1, i})
                    len = i-1;
                    break
                else
                    len = i;
                end
            end

            y = [];

            % nComp represents the number of outlet ports
            if opt.nZone == 4, nComp = 2; else nComp = 3; end
            for k = 1:nComp

                temp = cat(1, dyncData{k, 1:len});
                y = [y temp];

            end

            save('trajectory.dat', 'y' ,'-ascii', '-double', '-tabs');

        end % trajDataConvertToASCII

        function UIplot(str, opt, i, j, relativeDelta)
% ----------------------------------------------------------------------------------------
% Interactive plottings when enableDebug is enabled
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.UIplot: no enough input arguments \n');
            end

            if opt.enableDebug
                if strcmp('head', str)
                    fprintf(' ========================================================== \n');
                    fprintf(' ---- Interval ---- Switch ---- Cycle ---- CSS_relError ---- \n');
                elseif strcmp('interval', str)
                    if mod(j, opt.nInterval) == 1
                        fprintf(' %9d %12d \n', j, i);
                    else
                        fprintf(' %9d\n', j);
                    end
                elseif strcmp('cycle', str)
                    fprintf(' ---- Interval ---- Switch ---- Cycle ---- CSS_relError ---- \n');
                    fprintf(' %33d %17g \n', i/opt.nColumn, relativeDelta);
                end
            end

        end % UIplot

        function plotFigures(opt, plotData)
% ----------------------------------------------------------------------------------------
% This is the plot function
%
% Parameters:
% 		- opt. options
% 		- plotData. The data for plotting
% ----------------------------------------------------------------------------------------


            if nargin < 2
                error('SMB.plotFigures: There are no enough arguments \n');
            end

            if opt.enableDebug

                figure(01);clf

                y = plotData{1,1}.colState;
                for k = 2:opt.nColumn
                    y = [y; plotData{1,k}.colState];
                end

                FigSet = plot(y); axis([0, opt.nColumn*opt.nCellsColumn, 0,opt.yLim])
                ylabel('Concentration [mol/m^3]', 'FontSize', 10);
                if opt.nComponents == 2
                    legend('comp 1', 'comp 2', 'Location', 'NorthWest');
                elseif opt.nComponents == 3
                    legend('comp 1', 'comp 2', 'comp 3', 'Location', 'NorthWest');
                elseif opt.nComponents == 4
                    legend('comp 1', 'comp 2', 'comp 3', 'comp 4', 'Location', 'NorthWest');
                end

                set(FigSet, 'LineWidth', 2);
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                set(gca, 'ygrid', 'on');


                if opt.nZone == 4

                    if opt.nColumn == 4 && all( eq(opt.structID, ones(1,opt.nZone)) )

                        set(gca, 'XTick', (1/2:1:(opt.nColumn-0.5)).*opt.nCellsColumn);

                    elseif opt.nColumn == 8 && all( eq(opt.structID, ones(1,opt.nZone).*2) )

                        set(gca, 'XTick', (1:2:opt.nColumn).*opt.nCellsColumn);

                    elseif opt.nColumn == 12 && all( eq(opt.structID, ones(1,opt.nZone).*3) )

                        set(gca, 'XTick', (opt.nColumn/8 : 3: opt.nColumn).*opt.nCellsColumn);

                    elseif opt.nColumn == 16 && all( eq(opt.structID, ones(1,opt.nZone).*4) )

                        set(gca, 'XTick', (opt.nColumn/8 : 4: opt.nColumn).*opt.nCellsColumn);

                    end

                    set(gca, 'XTickLabel', {'Zone I','II','III','IV'});

                elseif opt.nZone == 5

                    if opt.nColumn == 5 && all( eq(opt.structID, ones(1,opt.nZone)) )

                        set(gca, 'XTick', (1/2:1:(opt.nColumn-0.5)).*opt.nCellsColumn);

                    elseif opt.nColumn == 10 && all( eq(opt.structID, ones(1,opt.nZone).*2) )

                        set(gca, 'XTick', (1:2:(opt.nColumn-1)).*opt.nCellsColumn);

                    elseif opt.nColumn == 15 && all( eq(opt.structID, ones(1,opt.nZone).*3) )

                        set(gca, 'XTick', (opt.nColumn/10 : 3: opt.nColumn).*opt.nCellsColumn);

                    elseif opt.nColumn == 20 && all( eq(opt.structID, ones(1,opt.nZone).*4) )

                        set(gca, 'XTick', (opt.nColumn/10 : 4: opt.nColumn).*opt.nCellsColumn);

                    end

                    set(gca, 'XTickLabel', {'Zone I','II','III','IV','V'});

                elseif opt.nZone == 8

                    if opt.nColumn == 8 && all( eq(opt.structID, ones(1,opt.nZone)) )

                        set(gca, 'XTick', (1/2:1:(opt.nColumn-0.5)).*opt.nCellsColumn);

                    elseif opt.nColumn == 16 && all( eq(opt.structID, ones(1,opt.nZone).*2) )

                        set(gca, 'XTick', (1:2:opt.nColumn).*opt.nCellsColumn);

                    elseif opt.nColumn == 24 && all( eq(opt.structID, ones(1,opt.nZone).*3) )

                        set(gca, 'XTick', (opt.nColumn/16:3:opt.nColumn).*opt.nCellsColumn);

                    elseif opt.nColumn == 32 && all( eq(opt.structID, ones(1,opt.nZone).*4) )

                        set(gca, 'XTick', (opt.nColumn/16:4:opt.nColumn).*opt.nCellsColumn);

                    end

                    set(gca, 'XTickLabel', {'Zone I','II','III','IV','V','VI','VII','VIII'});

                end % if opt.nZone

                for i = 1: (opt.nColumn-1)
                    line([i*opt.nCellsColumn,i*opt.nCellsColumn], [0, opt.yLim], 'color', 'k', 'LineStyle', '-.');
                end

                title(sprintf('%g / ', opt.structID));

            end % if opt.enableDebug

        end % plotFigures

        function plotDynamic(opt, dyncData, iter)
% ----------------------------------------------------------------------------------------
% This is the plot function
% The numbers in the figure() represent the number of the columns
% ----------------------------------------------------------------------------------------


            if nargin < 3
                error('SMB.plotDynamic: There are no enough arguments \n');
            end

            if opt.enableDebug

                figure(02);clf
                if opt.nZone == 4

                    for i = 1:2

                        y = cat(1, dyncData{i,:});

                        subplot(2,1,i);
                        FigSet = plot(y,'.'); axis([0, iter*opt.timePoints, 0,opt.yLim])
                        switch i
                            case 1
                                ylabel({'Raffinate Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                            case 2
                                ylabel({'Extract Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                        end
                        xString = sprintf('Switches in %3d-column case [n]', opt.nColumn);
                        xlabel(xString, 'FontSize', 10);

                        if opt.nComponents == 2 && i == 1
                            legend('comp 1', 'comp 2', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); set(FigSet(2),'Marker','*', 'MarkerSize',2);
                        elseif opt.nComponents == 3 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2);
                        elseif opt.nComponents == 4 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'comp 4', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','+', 'MarkerSize',2);
                        end

                        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                        set(gca, 'ygrid', 'on');
                        set(gca, 'XTick', opt.timePoints*(1:iter)*opt.nInterval);
                        set(gca, 'XTickLabel', (1:iter));


                        for j = 1: (iter-1)
                            line([j*opt.timePoints, j*opt.timePoints],[0,opt.yLim], 'color', 'k', 'LineStyle', ':');
                        end

                        for j = 1: (iter-1)
                            line([j*opt.timePoints*opt.nInterval, j*opt.timePoints*opt.nInterval],[0,opt.yLim], 'color', 'k', 'LineStyle', '-.');
                        end

                    end

                    suptitle('The dynamic trajectories at Raffinate and Extract ports');

                elseif opt.nZone == 5

                     for i = 1:3

                        y = cat(1, dyncData{i,:});

                        subplot(3,1,i);
                        FigSet = plot(y,'.'); axis([0, iter*opt.timePoints, 0,opt.yLim])
                        if strcmp('extract', opt.intermediate)
                            switch i
                                case 1
                                    ylabel({'Raffinate Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 2
                                    ylabel({'Extract_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 3
                                    ylabel({'Extract_1 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                            end
                        elseif strcmp('raffinate', opt.intermediate)
                            switch i
                                case 1
                                    ylabel({'Extract Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 2
                                    ylabel({'Raffinate_1 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 3
                                    ylabel({'Raffinate_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                            end
                        end
                        xString = sprintf('Switches in %3d-column case [n]', opt.nColumn);
                        xlabel(xString, 'FontSize', 10);

                        if opt.nComponents == 2 && i == 1
                            legend('comp 1', 'comp 2', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); set(FigSet(2),'Marker','*', 'MarkerSize',2);
                        elseif opt.nComponents == 3 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2);
                        elseif opt.nComponents == 4 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'comp 4', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','+', 'MarkerSize',2);
                        end

                        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                        set(gca, 'ygrid', 'on');
                        set(gca, 'XTick', opt.timePoints*(1:iter)*opt.nInterval);
                        set(gca, 'XTickLabel', (1:iter));


                        for j = 1: (iter-1)
                            line([j*opt.timePoints, j*opt.timePoints],[0,opt.yLim], 'color', 'k', 'LineStyle', ':');
                        end

                        for j = 1: (iter-1)
                            line([j*opt.timePoints*opt.nInterval, j*opt.timePoints*opt.nInterval],[0,opt.yLim], 'color', 'k', 'LineStyle', '-.');
                        end

                    end

                    suptitle('The dynamic trajectories at Raffinate and Extract ports');

                elseif opt.nZone == 8

                     for i = 1:3

                        y = cat(1, dyncData{i,:});

                        subplot(3,1,i);
                        FigSet = plot(y,'.'); axis([0, iter*opt.timePoints, 0,opt.yLim])
                        if strcmp('raffinate', opt.intermediate_feed)
                            switch i
                                case 1
                                    ylabel({'Raffinate_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 2
                                    ylabel({'Extract_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 3
                                    ylabel({'Extract_1 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                            end
                        elseif strcmp('extract', opt.intermediate_feed)
                            switch i
                                case 1
                                    ylabel({'Raffinate_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 2
                                    ylabel({'Extract_2 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                                case 3
                                    ylabel({'Raffinate_1 Port'; 'Concentration [mol/m^3]'}, 'FontSize', 10);
                            end
                        end
                        xString = sprintf('Switches in %3d-column case [n]', opt.nColumn);
                        xlabel(xString, 'FontSize', 10);

                        if opt.nComponents == 3 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2);
                        elseif opt.nComponents == 4 && i == 1
                            legend('comp 1', 'comp 2', 'comp 3', 'comp 4', 'Location', 'NorthWest');
                            set(FigSet(1),'Marker','^', 'MarkerSize',2); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',2); ...
                                set(FigSet(3),'Marker','+', 'MarkerSize',2);
                        end

                        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                        set(gca, 'ygrid', 'on');
                        set(gca, 'XTick', opt.timePoints*(1:iter)*opt.nInterval);
                        set(gca, 'XTickLabel', (1:iter));


                        for j = 1: (iter-1)
                            line([j*opt.timePoints, j*opt.timePoints],[0,opt.yLim], 'color', 'k', 'LineStyle', ':');
                        end

                        for j = 1: (iter-1)
                            line([j*opt.timePoints*opt.nInterval, j*opt.timePoints*opt.nInterval],[0,opt.yLim], 'color', 'k', 'LineStyle', '-.');
                        end

                    end

                    suptitle('The dynamic trajectories at Raffinate and Extract ports');

                end % if opt.nZone

            end % if opt.enableDebug

        end % plotDynamic


    end % methods


end % classdef
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
