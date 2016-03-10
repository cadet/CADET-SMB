
classdef SMB < handle
% =============================================================================
% This is the class of the functions of simulated moving bed.
% =============================================================================


    methods (Static = true, Access = 'public')


        function [outletProfile, lastState] = secColumn(inletProfile, params, lastState, ParSwarm)
% -----------------------------------------------------------------------------
% Simulation of the single column
%
% Parameters:
%       - inletProfile. Inlet time and corresponding concentration profile
%       - params. Get parameters for simulation
%       - lastState. The last STATE from previous simulation of next simulation
%       - ParSwarm. In the optimization situations, this is the vector which contains
%           the optimized decision variables
% 
% Returns:
%       - outletProfile. outlet time and corresponding concentration profile
%       - lastState. Record the last STATE which is used as the boundary condition
% -----------------------------------------------------------------------------


            if nargin < 4
                ParSwarm = [];
                if nargin < 3
                    lastState = [];
                    if nargin < 2
                        error('secColumn: There are no enough inputs for carrying out Simulator in CADET');
                    end
                end
            end

            if isempty(params.initMobilCon) && isempty(params.initSolidCon) && isempty(lastState)
                error('secColumn: There are no Initial / Boundary Conditions for the Simulator');
            end

%           Get parameters options
            [opt, ~, ~] = getParameters(ParSwarm);

            model = ModelGRM();
            model.nComponents = opt.nComponents;

%           The equilibrium isotherms
            if strcmp(opt.BindingModel, 'LinearBinding')

                model.kineticBindingModel = false;
                model.bindingModel = LinearBinding(); 

%               Adsorption parameters
                model.bindingParameters.LIN_KA   = opt.KA;
                model.bindingParameters.LIN_KD   = opt.KD;

            elseif strcmp(opt.BindingModel, 'MultiComponentLangmuirBinding')

                model.kineticBindingModel = true;
                model.bindingModel = MultiComponentLangmuirBinding();

                model.bindingParameters.MCL_KA   = opt.KA;
                model.bindingParameters.MCL_KD   = opt.KD;
                model.bindingParameters.MCL_QMAX = opt.QMAX;

            elseif strcmp(opt.BindingModel, 'MultiComponentBiLangmuirBinding')

                model.kineticBindingModel = true;
                model.bindingModel = MultiComponentBiLangmuirBinding();

                model.bindingParameters.MCL_KA1   = opt.KA(1);
                model.bindingParameters.MCL_KD1   = opt.KD(1);
                model.bindingParameters.MCL_QMAX1 = opt.QMAX(1);
                model.bindingParameters.MCL_KA2   = opt.KA(2);
                model.bindingParameters.MCL_KD2   = opt.KD(2);
                model.bindingParameters.MCL_QMAX2 = opt.QMAX(2);

            elseif strcmp(opt.BindingModel, 'StericMassAction')

                error('%s: it is not available yet.', opt.BindingModel);

            end

            if nargin >= 3 && ~isempty(lastState)
                model.initialState = lastState;
            else      
                model.initialMobileConcentration = params.initMobilCon;
                model.initialSolidConcentration  = params.initSolidCon;
            end

%           Transport
            model.dispersionColumn          = opt.dispersionColumn;
            model.filmDiffusion             = opt.filmDiffusion;
            model.diffusionParticle         = opt.diffusionParticle;
            model.diffusionParticleSurface  = opt.diffusionParticleSurface;
            model.interstitialVelocity      = params.interstitialVelocity;

%           Geometry
            model.columnLength        = opt.columnLength;
            model.particleRadius      = opt.particleRadius;
            model.porosityColumn      = opt.porosityColumn;
            model.porosityParticle    = opt.porosityParticle;

%           Apply the inlet profile to the CADET model
            Time = repmat({inletProfile.time}, 1, opt.nComponents);
            if opt.nComponents == 2
                Profile = [{inletProfile.concentration(:,1)}, {inletProfile.concentration(:,2)}];
            elseif opt.nComponents == 3
                Profile = [{inletProfile.concentration(:,1)}, {inletProfile.concentration(:,2)},...
                           {inletProfile.concentration(:,3)}];
            elseif opt.nComponents == 4
                Profile = [{inletProfile.concentration(:,1)}, {inletProfile.concentration(:,2)},...
                           {inletProfile.concentration(:,3)}, {inletProfile.concentration(:,4)}];
            end

            model.setInletsFromData(Time, Profile);

%           Turn off the warnings of the interpolation, optional
            warning('off', 'MATLAB:interp1:ppGriddedInterpolant');
            warning('off', 'MATLAB:interp1:UsePCHIP');

%           Discretization
            disc = DiscretizationGRM();
            disc.nCellsColumn   = opt.nCellsColumn;
            disc.nCellsParticle = opt.nCellsParticle;

%           Solving options
            sim = Simulator(model, disc);
            sim.nThreads = opt.nThreads;
            sim.solutionTimes = inletProfile.time;
            sim.solverOptions.time_integrator.ABSTOL         = opt.ABSTOL;
            sim.solverOptions.time_integrator.INIT_STEP_SIZE = opt.INIT_STEP_SIZE;
            sim.solverOptions.time_integrator.MAX_STEPS      = opt.MAX_STEPS;
            sim.solverOptions.WRITE_SOLUTION_ALL    = false;
            sim.solverOptions.WRITE_SOLUTION_LAST   = true;
            sim.solverOptions.WRITE_SENS_LAST       = false;
            sim.solverOptions.WRITE_SOLUTION_COLUMN_OUTLET = true;
            sim.solverOptions.WRITE_SOLUTION_COLUMN_INLET  = true;


%           Run the simulation
            try
                result = sim.simulate();
            catch e
                % Something went wrong
                error('CADET:simulationFailed', 'Check your settings and try again.\n%s',e.message);
            end

%           Extract the outlet profile
            outletProfile.time = result.solution.time;
            outletProfile.concentration = result.solution.outlet(:,:);
            lastState =  result.solution.lastState;


        end % secColumn

        function column = massConservation(currentData, interstVelocity, Feed, opt, sequence, alphabet)
% -----------------------------------------------------------------------------
% This is the function to calculate the concentration changes on each node.
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
% Fluid goes from Zone I to Zone II to Zone III, while the switch direction
% is from Zone I to Zone IV to Zone III;
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
% -----------------------------------------------------------------------------


            global stringSet;

%           Time points
            column.inlet.time = linspace(0, opt.switch, opt.timePoints);

%           Translate the alphabet into the position index of SMB unit
            index = SMB.stringIndexing(opt, alphabet);

            if ~strcmp(alphabet, 'a')
                pre_alphabet = char(alphabet - 1);
            else
                pre_alphabet = char(stringSet(opt.nColumn));
            end

%           Get the interstitial velocity of each column and boundary conditions
            params = SMB.getParams(sequence, interstVelocity, opt, index, alphabet, pre_alphabet);

            idx_i = eval( ['sequence' '.' alphabet] );     % the number of current column
            idx_j = eval( ['sequence' '.' pre_alphabet] ); % the number of the column before

%           Update the intersitial velocity, boundary conditions
            column.params = params{idx_i};

%           Calculate concentration of the column due to its position in the SMB unit
            switch index

                case 'D' % node DESORBENT

                    %   C_i^in = Q_{i-1} * C_{i-1}^out / Q_i
                    concentration = zeros(length(Feed.time), opt.nComponents);

                    column.inlet.concentration = concentration .* params{idx_j}.interstitialVelocity...
                        ./ params{idx_i}.interstitialVelocity; 

                case 'F' % node FEED

                    %   C_i^in = (Q_{i-1} * C_{i-1}^out + Q_F * C_F) / Q_i
                    column.inlet.concentration = (currentData{idx_j}.outlet.concentration .* ...
                        params{idx_j}.interstitialVelocity + Feed.concentration .* interstVelocity.feed) ...
                        ./ params{idx_i}.interstitialVelocity; 

                otherwise % node EXTRACT; RAFFINATE; MIDDLE

                    %   C_i^in = C_{i-1}^out
                    column.inlet.concentration = currentData{idx_j}.outlet.concentration;

            end

        end % massConservation

        function params = getParams(sequence, interstVelocity, opt, index, alphabet, pre_alphabet)
%-----------------------------------------------------------------------------------------
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
% 		- params. It contains interstitial velocity and boundary conditions
%-----------------------------------------------------------------------------------------


%             global stringSet;
%             string = char(stringSet(1:opt.nColumn));

            params = cell(1, opt.nColumn);
            for k = 1:opt.nColumn
%               set the initial conditions to the solver, but when lastState is used, this setup will be ignored 
                params{k} = struct('initMobilCon', zeros(1,opt.nComponents), 'initSolidCon',...
                    zeros(1,opt.nComponents), 'interstitialVelocity', []);
            end

%             for j = 1:opt.nColumn
% %               set the initial conditions to the solver, but when lastState is used, this setup will be ignored 
%                 params{ eval(['sequence' '.' string(j)]) }.initMobilCon = zeros(1,opt.nComponents);
%                 params{ eval(['sequence' '.' string(j)]) }.initSolidCon = zeros(1,opt.nComponents);
%             end

            if opt.nZone == 4

                switch index

                    case {'D' 'M_D'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                    case {'E' 'M_E'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle;
                    case {'F' 'M_F'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract;
                    case {'R' 'M_R'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract + interstVelocity.feed;

                end

            elseif opt.nZone == 5

                switch index

                    case {'D' 'M_D'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                    case {'E1' 'M_E1'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle;
                    case {'E2' 'M_E2'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1 - interstVelocity.extract2;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1;
                    case {'F' 'M_F'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent + interstVelocity.raffinate;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.extract1 - interstVelocity.extract2;
                    case {'R' 'M_R'}
                        params{eval(['sequence' '.' alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent;
                        params{eval(['sequence' '.' pre_alphabet])}.interstitialVelocity = interstVelocity.recycle - interstVelocity.desorbent + interstVelocity.raffinate;

                end

            end

        end % getParams


        function index = stringIndexing(opt, alphabet)
%-----------------------------------------------------------------------------------------
% This is the function which interpret the alphabet of columns into the position of
% SMB unit.
% 
% Parameters:
% 		- opt. options involving the parameters for models
% 		- alphabet. The letter is used for switching
%
% Returns:
% 		- index. The indexing letter used for calculation of the mass conservation
%-----------------------------------------------------------------------------------------


            global stringSet;

            string = stringSet(1:opt.nColumn);

%           Preallocation
%           stringBlock is used for storing the alphabet in each zone
            stringBlock = cell(1, opt.nZone);

%           Separate the string into nZone cells
            stringBlock{1} = string(1:opt.structID(1));
            for k = 2:opt.nZone
                stringBlock{k} = string( sum(opt.structID(1:k-1))+1 : sum(opt.structID(1:k)) );
            end

%           Assign each alphabet with the indexing letter, D,E,F,R
            if opt.nZone == 4

                if any( strcmp(alphabet, stringBlock{1}) )
                    if strcmp(alphabet, stringBlock{1}(1))
                        index = 'D';
                    else
                        index = 'M_D';
                    end

                elseif any( strcmp(alphabet, stringBlock{2}) )
                    if strcmp(alphabet, stringBlock{2}(1))
                        index = 'E';
                    else
                        index = 'M_E';
                    end

                elseif any( strcmp(alphabet, stringBlock{3}) )
                    if strcmp(alphabet, stringBlock{3}(1))
                        index = 'F';
                    else
                        index = 'M_F';
                    end

                elseif any( strcmp(alphabet, stringBlock{4}) )
                    if strcmp(alphabet, stringBlock{4}(1))
                        index = 'R';
                    else
                        index = 'M_R';
                    end

                end

            elseif opt.nZone == 5

                if any( strcmp(alphabet, stringBlock{1}) )
                    if strcmp(alphabet, stringBlock{1}(1))
                        index = 'D';
                    else
                        index = 'M_D';
                    end

                elseif any( strcmp(alphabet, stringBlock{2}) )
                    if strcmp(alphabet, stringBlock{2}(1))
                        index = 'E1';
                    else
                        index = 'M_E1';
                    end

                elseif any( strcmp(alphabet, stringBlock{3}) )
                    if strcmp(alphabet, stringBlock{3}(1))
                        index = 'E2';
                    else
                        index = 'M_E2';
                    end

                elseif any( strcmp(alphabet, stringBlock{4}) )
                    if strcmp(alphabet, stringBlock{4}(1))
                        index = 'F';
                    else
                        index = 'M_F';
                    end

                elseif any( strcmp(alphabet, stringBlock{5}) )
                    if strcmp(alphabet, stringBlock{5}(1))
                        index = 'R';
                    else
                        index = 'M_R';
                    end

                end

            end

        end % stringIndexing

        function obj = positionIndexing(opt)
%-----------------------------------------------------------------------------------------
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
%-----------------------------------------------------------------------------------------


            Number = (1:opt.nColumn);

%           numberBlock is used for storing the column number in each zone
            numberBlock = cell(1, opt.nZone);

%           separate the number string into nZone cells
            numberBlock{1} = Number(1:opt.structID(1));
            for k = 2:opt.nZone
                numberBlock{k} = Number( sum(opt.structID(1:k-1))+1 : sum(opt.structID(1:k)) );
            end

            if opt.nZone == 4

                obj.position_ext = numberBlock{1}(end); % Desorbent node
                obj.position_raf = numberBlock{3}(end); % Feed node

            elseif opt.nZone == 5

                obj.position_ext1 = numberBlock{1}(end); % Desorbent node
                obj.position_ext2 = numberBlock{2}(end); % Extract_1 node
                obj.position_raf = numberBlock{4}(end);  % Feed node

            end

        end % positionIndexing


        function objective = objectiveFunction(Results, opt)
%-----------------------------------------------------------------------------------------
% The objective function for the optimizers
% You can also define your own objective function here. The default function is: 
%
% Max Productivity_extract + Productivity_raffinate
% s.t. Purity_extract   >= 99% for more retained component
%      Purity_raffinate >= 99% for less retained component
%      other implicit constraints, such as upbound on Desorbent consumption
%-----------------------------------------------------------------------------------------


            if opt.nZone == 4
%               Construct the Penalty Function for the objective function
                penalty = abs( min(Results.Purity_extract - opt.Purity_extract_limit, 0) ) * opt.Penalty_factor ...
                    + abs( min(Results.Purity_raffinate - opt.Purity_raffinate_limit, 0) ) * opt.Penalty_factor;

%               (-) since in the optimizer, the defined program is of optimization of minimum.    
                objective = -(Results.Productivity_extract + Results.Productivity_raffinate) + penalty;

            elseif opt.nZone == 5
%               Construct the Penalty Function for the objective function
                penalty = abs( min(Results.Purity_extract1 - opt.Purity_extract1_limit, 0) ) * opt.Penalty_factor ...
                    + abs( min(Results.Purity_extract2 - opt.Purity_extract2_limit, 0) ) * opt.Penalty_factor ...
                    + abs( min(Results.Purity_raffinate - opt.Purity_raffinate_limit, 0) ) * opt.Penalty_factor;

%               (-) since in the optimizer, the defined program is of optimization of minimum.    
                objective = -(Results.Productivity_extract1 + Results.Productivity_extract2 + Results.Productivity_raffinate) + penalty;
            end


            if opt.enableDebug
                fprintf('**** The objective value:  %g \n', objective);
            end

        end % objectiveFunction

        function Results = Purity_Productivity(plotData, opt)
%-----------------------------------------------------------------------------------------
% Calculation of the performance index of SMB, such Purity and Productivity
%
% Parameters:
% 		- plotData. The data for plotting
% 		- opt. options of parameters
%
% Returns:
% 		- Resutls. A struct data which contains the purity and productivity of extract
% 			and raffinate ports, respectively
%-----------------------------------------------------------------------------------------


%           Get the position of the withdrawn ports to calculate the purity and productivity
            obj = SMB.positionIndexing(opt);

%           The nominator of the formualar of productivity
            Nominator = pi * (opt.columnDiameter/2)^2 * opt.columnLength * (1-opt.porosityColumn);

            if opt.nZone == 4

%               using column 1 to calculate the integral of purity, plotData{1,x}
                position_ext = obj.position_ext; position_raf = obj.position_raf;

%               Please be quite careful, which component is used for statistics (change them with comp_ext_ID or comp_raf_ID)
                sum_ext = 0; sum_raf = 0;
                for k = 1:opt.nComponents
                    sum_ext = sum_ext + trapz(plotData{1,position_ext}.outlet.time, plotData{1,position_ext}.outlet.concentration(:,k));
                    sum_raf = sum_raf + trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,k));
                end

%               Extract ports
                Purity_extract = trapz(plotData{1,position_ext}.outlet.time, plotData{1,position_ext}.outlet.concentration(:,opt.comp_ext_ID)) / sum_ext;

%               Raffinate ports
                Purity_raffinate = trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,opt.comp_raf_ID)) / sum_raf;

%               per switching time, in the tank of extract port, such (unit: g/m^3) amount of target component was collected.
                Productivity_extract = trapz(plotData{1,position_ext}.outlet.time, plotData{1,position_ext}.outlet.concentration(:,opt.comp_ext_ID))...
                    * opt.molMass(opt.comp_ext_ID) * opt.flowRate_extract / Nominator;

                Productivity_raffinate = trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,opt.comp_raf_ID))...
                    * opt.molMass(opt.comp_raf_ID) * opt.flowRate_raffinate / Nominator;


                if opt.enableDebug
                    fprintf('Purity (Extract): %g %% \n', Purity_extract * 100);
                    fprintf('Purity (Raffinate): %g %% \n', Purity_raffinate * 100)
                    fprintf('Productivity (Extract) in each switching time: %g g/m^3 \n', Productivity_extract);
                    fprintf('Productivity (Raffinate) in each switching time: %g g/m^3 \n', Productivity_raffinate);
                end

                Results = struct('Purity_extract', Purity_extract, 'Purity_raffinate', Purity_raffinate, ...
                    'Productivity_extract', Productivity_extract, 'Productivity_raffinate', Productivity_raffinate);
                Results.Data = plotData;


            elseif opt.nZone == 5

                position_ext1 = obj.position_ext1; position_ext2 = obj.position_ext2; position_raf = obj.position_raf;

%               Please be quite careful, which component is used for statistics (change them with comp_ext_ID or comp_raf_ID)
                sum_ext1 = 0; sum_ext2 = 0; sum_raf = 0;
                for k = 1:opt.nComponents
                    sum_ext1 = sum_ext1 + trapz(plotData{1,position_ext1}.outlet.time, plotData{1,position_ext1}.outlet.concentration(:,k));
                    sum_ext2 = sum_ext2 + trapz(plotData{1,position_ext2}.outlet.time, plotData{1,position_ext2}.outlet.concentration(:,k));
                    sum_raf = sum_raf + trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,k));
                end

%               Extract ports
                Purity_extract1 = trapz(plotData{1,position_ext1}.outlet.time, plotData{1,position_ext1}.outlet.concentration(:,opt.comp_ext1_ID)) / sum_ext1;
                Purity_extract2 = trapz(plotData{1,position_ext2}.outlet.time, plotData{1,position_ext2}.outlet.concentration(:,opt.comp_ext2_ID)) / sum_ext2;

%               Raffinate ports
                Purity_raffinate = trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,opt.comp_raf_ID)) / sum_raf;

%               per switching time, in the tank of extract port, such (unit: g/m^3) amount of target component was collected.
                Productivity_extract1 = trapz(plotData{1,position_ext1}.outlet.time, plotData{1,position_ext1}.outlet.concentration(:,opt.comp_ext1_ID))...
                    * opt.molMass(opt.comp_ext1_ID) * opt.flowRate_extract1 / Nominator;

                Productivity_extract2 = trapz(plotData{1,position_ext2}.outlet.time, plotData{1,position_ext2}.outlet.concentration(:,opt.comp_ext2_ID))...
                    * opt.molMass(opt.comp_ext2_ID) * opt.flowRate_extract2 / Nominator;

                Productivity_raffinate = trapz(plotData{1,position_raf}.outlet.time, plotData{1,position_raf}.outlet.concentration(:,opt.comp_raf_ID))...
                    * opt.molMass(opt.comp_raf_ID) * opt.flowRate_raffinate / Nominator;


                if opt.enableDebug
                    fprintf('Purity (Extract_1): %g %% \n', Purity_extract1 * 100);
                    fprintf('Purity (Extract_2): %g %% \n', Purity_extract2 * 100);
                    fprintf('Purity (Raffinate): %g %% \n', Purity_raffinate * 100)
                    fprintf('Productivity (Extract_1) in each switching time: %g g/m^3 \n', Productivity_extract1);
                    fprintf('Productivity (Extract_2) in each switching time: %g g/m^3 \n', Productivity_extract2);
                    fprintf('Productivity (Raffinate) in each switching time: %g g/m^3 \n', Productivity_raffinate);
                end

                Results = struct('Purity_extract1', Purity_extract1, 'Purity_extract2', Purity_extract2,...
                    'Purity_raffinate', Purity_raffinate, 'Productivity_extract1', Productivity_extract1,...
                    'Productivity_extract2', Productivity_extract2, 'Productivity_raffinate', Productivity_raffinate);
                Results.Data = plotData;

            end

        end % Purity_Productivity


        function columnSpline = CSTR(Profile, column, opt)
% -----------------------------------------------------------------------------
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
% -----------------------------------------------------------------------------


            for i = 1:opt.nComponents

                IC = Profile.concentration(1,i);

                [T, Y] = ode15s(@(t, y) ode_CSTR(t, y, Profile.time, Profile.concentration, column.params.interstitialVelocity, i), [0, opt.switch], IC);

                columnSpline.concentration(:,i) = interp1(T, Y, Profile.time);

            end

            columnSpline.time = Profile.time;


            function DyDt = ode_CSTR(t, y, time, concentration, interstVelocity, comp)
            % -----------------------------------------------------------------
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
            % -----------------------------------------------------------------


                tankLength = 0.01;

                concentration_f = interp1(time, concentration(:,comp), t);

                tau = tankLength / interstVelocity;

                DyDt = 1/tau * (concentration_f - y);

            end

        end % CSTR

        function [columnSpline, lastState]  = DPFR(Profile, initialState, opt)
% -----------------------------------------------------------------------------
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
% -----------------------------------------------------------------------------


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
            % -----------------------------------------------------------------
            % The PDE function of dispersive plug flow reactor
            % Dy/Dt = - u Dy/Dz + D_{ax} DDy/DDz
            %
            % Parameters:
            %       - t. Time point
            %       - x. Initial conditons
            % -----------------------------------------------------------------


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


        function plotFigures(opt, plotData)
%-----------------------------------------------------------------------------------------
% This is the plot function 
% The numbers in the figure() represent the number of the columns
% 
% Parameters:
% 		- opt. options
% 		- plotData. The data for plotting
%-----------------------------------------------------------------------------------------


            if nargin < 2
                disp('Error in plotFigures: there are no enough input data for the function');
            else
                if isempty(opt)
                    disp('Error in plotFigures: the options of the parameters are missing');
                elseif isempty(plotData)
                    disp('Error in plotFigures: the data for figures plotting are missing');
                end
            end

            if opt.enableDebug

                figure(01);clf

                y = zeros(1, opt.nComponents);
                for k = opt.nColumn:-1:1
                    y = [y; plotData{k}.outlet.concentration];
                end

                FigSet = plot(y); axis([0,opt.nColumn*opt.timePoints, 0,opt.yLim])
                ylabel('Concentration [Mol]', 'FontSize', 10);
                if opt.nComponents == 2
                    legend('comp 1', 'comp 2');
                elseif opt.nComponents == 3
                    legend('comp 1', 'comp 2', 'comp 3');
                elseif opt.nComponents == 4
                    legend('comp 1', 'comp 2', 'comp 3', 'comp 4');
                end

                set(FigSet, 'LineWidth', 2);
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                set(gca, 'ygrid', 'on');


                if opt.nZone == 4

                    if opt.nColumn == 4 && all( eq(opt.structID, [1 1 1 1]) )

                        set(gca, 'XTick', (1/2:1:(opt.nColumn-0.5)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 8 && all( eq(opt.structID, [2 2 2 2]) )

                        set(gca, 'XTick', (1:2:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 12 && all( eq(opt.structID, [3 3 3 3]) )

                        set(gca, 'XTick', (opt.nColumn/8:3:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 16 && all( eq(opt.structID, [4 4 4 4]) )

                        set(gca, 'XTick', (opt.nColumn/8:4:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone IV','Zone III','Zone II','Zone I'});

                    end

                elseif opt.nZone == 5

                    if opt.nColumn == 5 && all( eq(opt.structID, [1 1 1 1 1]) )

                        set(gca, 'XTick', (1/2:1:(opt.nColumn-0.5)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone V','Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 10 && all( eq(opt.structID, [2 2 2 2 2]) )

                        set(gca, 'XTick', (1:2:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone V','Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 15 && all( eq(opt.structID, [3 3 3 3 3]) )

                        set(gca, 'XTick', (opt.nColumn/10:3:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone V','Zone IV','Zone III','Zone II','Zone I'});

                    elseif opt.nColumn == 20 && all( eq(opt.structID, [4 4 4 4 4]) )

                        set(gca, 'XTick', (opt.nColumn/10:4:(opt.nColumn-1)).*opt.timePoints);
                        set(gca, 'XTickLabel', {'Zone V','Zone IV','Zone III','Zone II','Zone I'});

                    end

                end % if opt.nZone == 4 / opt.nZone == 5

                for i = 1: (opt.nColumn-1)
                    line([i*opt.timePoints,i*opt.timePoints], [0, opt.yLim], 'color', 'k', 'LineStyle', '-.');
                end

            end % if opt.enableDebug

        end % plotFigures

        function plotDynamic(opt, dyncData, iter)
%-----------------------------------------------------------------------------------------
%  This is the plot function
%  The numbers in the figure() represent the number of the columns
%-----------------------------------------------------------------------------------------


            if nargin < 2
                disp('Error: there are no enough input data for the function, plotFigures');
            else
                if isempty(opt)
                    disp('Error in plotFigures: the options of the parameters are missing');
                elseif isempty(dyncData)
                    disp('Error in plotFigures: the data for figures plotting are missing');
                end
            end

            if opt.enableDebug

                figure(100);clf
                if opt.nZone == 4

                    for i = 1:2

                        y = cat(1, dyncData{i,:});

                        subplot(2,1,i);
                        FigSet = plot(y,'.'); axis([0,iter*opt.timePoints, 0,opt.yLim])
                        ylabel('Concentration [Mol]', 'FontSize', 10);

                        if opt.nComponents == 2
                            legend('comp 1', 'comp 2');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); set(FigSet(2),'Marker','*', 'MarkerSize',3.5);
                        elseif opt.nComponents == 3
                            legend('comp 1', 'comp 2', 'comp 3');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',3.5);
                        elseif opt.nComponents == 4
                            legend('comp 1', 'comp 2', 'comp 3', 'comp 4');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','+', 'MarkerSize',3.5);
                        end

                        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                        set(gca, 'ygrid', 'on');
                        set(gca, 'XTick', size(y,1)/2);

                        switch i
                            case 1
                                set(gca, 'XTickLabel', {'Raffinate Port'});
                            case 2
                                set(gca, 'XTickLabel', {'Extract Port'});
                        end

                        for j = 1: (iter-1)
                            line([j*opt.timePoints, j*opt.timePoints],[0,opt.yLim], 'color', 'k', 'LineStyle', '-.');
                        end

                    end

                    suptitle('The evolution of the concentration from the Raffinate port and Extract port');

                elseif opt.nZone == 5

                     for i = 1:3

                        y = cat(1, dyncData{i,:});

                        subplot(3,1,i);
                        FigSet = plot(y,'.'); axis([0,iter*opt.timePoints, 0,opt.yLim])
                        ylabel('Concentration [Mol]', 'FontSize', 10);

                        if opt.nComponents == 2
                            legend('comp 1', 'comp 2');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); set(FigSet(2),'Marker','*', 'MarkerSize',3.5);
                        elseif opt.nComponents == 3
                            legend('comp 1', 'comp 2', 'comp 3');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',3.5);
                        elseif opt.nComponents == 4
                            legend('comp 1', 'comp 2', 'comp 3', 'comp 4');
                            set(FigSet(1),'Marker','^', 'MarkerSize',3.5); ...
                                set(FigSet(2),'Marker','*', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','s', 'MarkerSize',3.5); ...
                                set(FigSet(3),'Marker','+', 'MarkerSize',3.5);
                        end

                        set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
                        set(gca, 'ygrid', 'on');
                        set(gca, 'XTick', size(y,1)/2);

                        switch i
                            case 1
                                set(gca, 'XTickLabel', {'Raffinate Port'});
                            case 2
                                set(gca, 'XTickLabel', {'Extract_2 Port'});
                            case 3
                                set(gca, 'XTickLabel', {'Extract_1 Port'});
                        end

                        for j = 1: (iter-1)
                            line([j*opt.timePoints, j*opt.timePoints],[0,opt.yLim], 'color', 'k', 'LineStyle', '-.');
                        end

                    end

                    suptitle('The evolution of the concentration from the Raffinate port and Extract ports');

                end

            end

        end % plotDynamic


    end % methods


end % classdef
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