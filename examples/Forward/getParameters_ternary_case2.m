function [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin)
% =============================================================================
%   Case 2, a 5-zone ten-column case for ternary separation (two extract)
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


    % The parameter setting for simulator
    opt.tolIter         = 1e-4;  % tolerance of the SMB stopping criterion
    opt.nMaxIter        = 1000;  % the maximum iteration step in SMB
    opt.nThreads        = 4;     % threads of CPU, up to your computer
    opt.nCellsColumn    = 40;    % discretization number in one column
    opt.nCellsParticle  = 1;     % discretization number in one particle
    opt.ABSTOL          = 1e-10; % tolerance of CADET stopping criterion
    opt.INIT_STEP_SIZE  = 1e-14; % refer your to CADET manual
    opt.MAX_STEPS       = 5e6;   % the maximum iteration step in CADET
    opt.enableDebug     = true;  % set it true when you want to see the figures

    % The parameter setting for the SMB
    opt.switch          = 1394.4;   % switching time (s)
    opt.timePoints      = 1000;  % the observed time-points
    opt.Purity_limit    = [0.95, 0.50, 0.99];  % used for constructing constraints
    opt.Penalty_factor  = 10;    % penalty factor in penalty function

    % Network configuration
    opt.nZone       = 5;
    opt.nColumn     = 10;
    opt.structID    = [2 2 2 2 2]; % 10 columns in five-zone, 2 in each zone
    opt.intermediate = 'extract'; % two extract configuration

    % Binding: Linear Binding isotherm
    opt.BindingModel = 'LinearBinding';
    opt.nComponents = 3;
    opt.KA = [5.34, 6.80, 11.20]; % [comp_A, comp_B, comp_C], A,B for raffinate, C for extract
    opt.KD = [1, 1, 1];   % K_A < K_B < K_C
    opt.compTargID  = [3, 2, 1]; % target components at [Extract1 Extract2 Raffinate] ports

    % Transport
    opt.dispersionColumn          = ones(1, opt.nZone) .* 1.1781e-7;    % D_{ax} m^2/s
    opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];           % K_f m/s
    opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];              % D_p m^2/s
    opt.diffusionParticleSurface  = [0.0, 0.0, 0.0];

    % Geometry
    opt.columnLength        = 25e-2;     % m
    opt.columnDiameter      = 1.0e-2;    % m
    opt.particleRadius      = 15e-6/2;   % m
    opt.porosityColumn      = 0.47;
    opt.porosityParticle    = 0.0000001; % e_p very small to ensure e_t = e_c

    % Parameter units transformation
    % The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;  % m^2
    flowRate.recycle    = 9.7333e-8;            % m^3/s
    flowRate.feed       = 4.1667e-9;            % m^3/s
    flowRate.raffinate  = 1.2000e-8;            % m^3/s
    flowRate.desorbent  = 5.6993e-8;            % m^3/s
    flowRate.extract1   = 3.7833e-8;            % m^3/s
    flowRate.extract2   = 1.1333e-8;            % m^3/s
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

    concentrationFeed   = [1.0, 1.0, 1.0];    % g/m^3
    opt.molMass         = [309.401, 309.401, 309.401]; % g/mol
    opt.yLim            = max(concentrationFeed ./ opt.molMass) * 1.1; % mol/m^3

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
