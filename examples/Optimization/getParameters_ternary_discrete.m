function [opt, interstVelocity, Feed, Desorbent] = getParameters(ParSwarm)
% =============================================================================
%   Case #2, a 5-zone case for column configuration optimization
%
% This is the function to input all the necessary data for simulation
% Returns:
%       1. opt stands for options, which involves the parameter settings
%       for the algorithm, the binding isotherm, and the model equations
%
%       2. interstVelocity is calculated from flowrate of each column and inlet.
%       interstitial_velocity = flow_rate / (across_area * porosity_Column)
%
%       3. Feed initializes the injection concentration
% =============================================================================


    valueAssign = struct('columnLength',ParSwarm(1), 'switch',ParSwarm(2),...
        'recycle',ParSwarm(3),'feed',ParSwarm(4), 'desorbent',ParSwarm(5),...
        'extract1',ParSwarm(6),'extract2',ParSwarm(7));

    % The parameter setting for simulator
    opt.tolIter         = 1e-3;   % tolerance of the SMB stopping criterion
    opt.nMaxIter        = 1000;   % the maximum iteration step in SMB
    opt.nThreads        = 4;      % threads of CPU, up to your computer
    opt.nCellsColumn    = 30;     % discretization number in one column
    opt.nCellsParticle  = 1;      % discretization number in one particle
    opt.ABSTOL          = 1e-9;   % tolerance of CADET stopping criterion
    opt.INIT_STEP_SIZE  = 1e-14;  % refer your to CADET manual
    opt.MAX_STEPS       = 5e6;    % the maximum iteration step in CADET
    opt.enableDebug     = false;  % set it false if you are using the optimizer

    % The parameter setting for the SMB
    opt.switch          = valueAssign.switch;  % switching time (s)
    opt.timePoints      = 1000;         % the observed time-points
    opt.Purity_limit    = [0.9750, 0.6946, 0.9995];  % used for constructing constraints
    opt.Penalty_factor  = 10;    % penalty factor in penalty function

    % Network configuration
    opt.nZone   = 5;
    opt.nColumn = 7; % 7 columns in five-zone
    opt.structNumber = 20;
    opt.intermediate = 'extract'; % two extract configuration

    % Binding: Linear Binding isotherm
    opt.BindingModel = 'LinearBinding';
    opt.nComponents = 3;
    opt.KA = [3.15, 7.4, 23]; % [comp_A, comp_B], A for raffinate, B for extract
    opt.KD = [1, 1, 1];
    opt.compTargID  = [3, 2, 1]; % target components at [Extract1 Extract2 Raffinate] ports

    % Transport
    opt.dispersionColumn          = ones(1, opt.nZone) .* 3.8148e-10;   % D_{ax} m^2/s
    opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];           % K_f m/s
    opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];              % D_p m^2/s
    opt.diffusionParticleSurface  = [0.0, 0.0, 0.0];

    % Geometry
    opt.columnLength        = valueAssign.columnLength; % m
    opt.columnDiameter      = 1.0e-2;     % m
    opt.particleRadius      = 30e-6/2;    % m
    opt.porosityColumn      = 0.8;
    opt.porosityParticle    = 0.00000001; % e_p very small to ensure e_t = e_c

    % Parameter units transformation
    % The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;      % m^2
    flowRate.recycle    = valueAssign.recycle;      % m^3/s
    flowRate.feed       = valueAssign.feed;         % m^3/s
    flowRate.desorbent  = valueAssign.desorbent;    % m^3/s
    flowRate.extract1   = valueAssign.extract1;     % m^3/s
    flowRate.extract2   = valueAssign.extract2;     % m^3/s
    flowRate.raffinate  = flowRate.desorbent - flowRate.extract1 - flowRate.extract2 + flowRate.feed;  % m^3/s
    opt.flowRate_extract1  = flowRate.extract1;
    opt.flowRate_extract2  = flowRate.extract2;
    opt.flowRate_raffinate = flowRate.raffinate;

    % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract1  = flowRate.extract1 / (crossArea*opt.porosityColumn);     % m/s
    interstVelocity.extract2  = flowRate.extract2 / (crossArea*opt.porosityColumn);     % m/s

    concentrationFeed   = [1.0, 1.0, 1.0];   % g/m^3
    opt.molMass         = [227.217, 267.24, 251.24192]; % g/mol
    opt.yLim            = max(concentrationFeed ./ opt.molMass); % mol/m^3

    % Feed concentration setup
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), opt.nComponents);
    Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);

    for i = 1:opt.nComponents
        Feed.concentration(1:end, i) = concentrationFeed(i) / opt.molMass(i); % mol/m^3
    end

% -----------------------------------------------------------------------------
    % Capable of placing a CSTR or DPFR apparatues before and after the calculated column

    % Continuous Stirred Tank Reactor
    opt.enable_CSTR = false;
    opt.CSTR_length = 0.01;

    % Dispersive Plug Flow Reactor
    opt.enable_DPFR = false;

    opt.DPFR_length = 0.0066;
    opt.DPFR_nCells = 50;

    opt.DPFR_velocity   = 0.00315;
    opt.DPFR_dispersion = 2.5e-20;


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
