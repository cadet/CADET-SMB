function [opt, interstVelocity, Feed] = getParameters(ParSwarm)
%   Case 5, a five-column demonstration case for ternary separation

% =============================================================================
% This is the function to input all the necessary data for simulation
%
% Returns: 
%       1. opt stands for options, which involves the parameter settings
%       for the algorithm, the binding isotherm, and the model equations
%
%       2. interstVelocity is calculated from flowrate of each column and inlet. 
%       interstitial_velocity = flow_rate / (across_area * porosity_Column)
%
%       3. Feed initializes the injection concentration
% =============================================================================


    valueAssign = struct('columnLength',ParSwarm(1), 'switch',ParSwarm(2), 'recycle',ParSwarm(3),...
        'feed',ParSwarm(4), 'desorbent',ParSwarm(5), 'extract',ParSwarm(6));

%   The parameter setting for simulator
    opt.tolIter         = 1e-3;   % tolerance of the SMB stopping criterion
    opt.nMaxIter        = 1000;   % the maximum iteration step in SMB
    opt.nThreads        = 8;      % threads of CPU, up to your computer
    opt.nCellsColumn    = 30;     % discretization number in one column
    opt.nCellsParticle  = 1;      % discretization number in one particle
    opt.ABSTOL          = 1e-9;   % tolerance of CADET stopping criterion
    opt.INIT_STEP_SIZE  = 1e-14;  % refer your to CADET manual
    opt.MAX_STEPS       = 5e6;    % the maximum iteration step in CADET

%   The parameter setting for the SMB
    opt.switch          = valueAssign.switch;  % s % switching time 
    opt.timePoints      = 1000;         % the observed time-points
    opt.Purity_extract1_limit  = 0.95;  % used for constructing constraints
    opt.Purity_extract2_limit  = 0.69;  % used for constructing constraints
    opt.Purity_raffinate_limit = 0.99;  % used for constructing constraints
    opt.Penalty_factor         = 10;    % penalty factor in penalty function

    opt.enableDebug = false; % set it false if you are using the optimizer
    opt.nZone   = 5;
    opt.nColumn = 5;         % 4,8,12,16 -column cases are available

%   Binding: Linear Binding isotherm
    opt.BindingModel = 'LinearBinding';
    opt.nComponents = 3;
    opt.KA = [3.15, 7.4, 23]; % [comp_A, comp_B], A for raffinate, B for extract
    opt.KD = [1, 1, 1];
    opt.comp_raf_ID  = 1; % the target component withdrawn from the raffinate ports
    opt.comp_ext1_ID = 3; % the target component withdrawn from the extract ports
    opt.comp_ext2_ID = 2; % the target component withdrawn from the extract ports

%   Transport
    opt.dispersionColumn          = 3.8148e-10;     % D_{ax}
    opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];      % K_f 
    opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];  % D_p
    opt.diffusionParticleSurface  = [0.0, 0.0, 0.0];

%   Geometry
    opt.columnLength        = valueAssign.columnLength;      % m
    opt.columnDiameter      = 1.0e-2;     % m
    opt.particleRadius      = 30e-6/2;    % m % user-defined one in this case
    opt.porosityColumn      = 0.8;
    opt.porosityParticle    = 0.00000001;  % unknown

%   Parameter units transformation
%   The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;      % m^2
    flowRate.recycle    = valueAssign.recycle;      % m^3/s  
    flowRate.feed       = valueAssign.feed;         % m^3/s
    flowRate.desorbent  = valueAssign.desorbent;    % m^3/s
    flowRate.extract1   = valueAssign.extract;      % m^3/s
    flowRate.extract2   = 4.6367e-8;                % m^3/s
    flowRate.raffinate  = flowRate.desorbent - flowRate.extract1 - flowRate.extract2 + flowRate.feed;  % m^3/s
    opt.flowRate_extract1  = flowRate.extract1;
    opt.flowRate_extract2  = flowRate.extract2;
    opt.flowRate_raffinate = flowRate.raffinate;

%   Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s 
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract1  = flowRate.extract1 / (crossArea*opt.porosityColumn);      % m/s
    interstVelocity.extract2  = flowRate.extract2 / (crossArea*opt.porosityColumn);      % m/s

    concentrationFeed 	= [1.0, 1.0, 1.0];   % g/m^3 [concentration_compA, concentration_compB]
    opt.molMass         = [227.217, 267.24, 251.24192]; % The molar mass of each components
    opt.yLim            = max(concentrationFeed ./ opt.molMass); % the magnitude for plotting

%   Feed concentration setup   
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), opt.nComponents);

    for i = 1:opt.nComponents
        Feed.concentration(1:end,i) = (concentrationFeed(i) / opt.molMass(i));
    end

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, either binary or ternary.
%  
%  Author: QiaoLe He   E-mail: q.he@fz-juelich.de
%                                      
%  Institute: Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. Please see the license of CADET.
% =============================================================================