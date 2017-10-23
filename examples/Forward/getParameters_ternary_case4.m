function [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin)
% =============================================================================
%   Case 4, a 8-zone eight-column case for ternary separation
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
    opt.tolIter         = 1e-4;
    opt.nMaxIter        = 1000;
    opt.nThreads        = 4;
    opt.nCellsColumn    = 40;
    opt.nCellsParticle  = 1;
    opt.ABSTOL          = 1e-10;
    opt.INIT_STEP_SIZE  = 1e-14;
    opt.MAX_STEPS       = 5e6;

    % The parameter setting for the SMB
	opt.nInterval 		= 20;
    opt.switch          = 1552/opt.nInterval;
    opt.timePoints      = 1000/opt.nInterval;
    opt.Purity_limit    = [0.50, 0.95, 0.50];
    opt.Penalty_factor  = 10;

    opt.enableDebug = true;
    opt.nZone       = 8;    % 8-zone for ternary separation
    opt.nColumn     = 8;
    opt.structID    = [1 1 1 1 1 1 1 1];
    opt.intermediate_feed = 'raffinate';

    % Binding: Linear Binding isotherm
    opt.BindingModel = 'LinearBinding';
    opt.nComponents = 3;
    opt.KA = [0.26 0.43 0.6]; % [comp_A, comp_B, comp_C], A for raffinate2, B for extract2, C for extract1
    opt.KD = [1, 1, 1];
    opt.compTargID = [3, 2, 1]; % target components at [Extract1 Extract2 Raffinate] ports

    % Transport
    opt.dispersionColumn          = ones(1, opt.nZone) .* 3.8148e-16; % D_{ax}
    opt.filmDiffusion             = [5.0e-5, 5.0e-5, 5.0e-5]; % K_f
    opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];    % D_p
    opt.diffusionParticleSurface  = [0.0, 0.0, 0.0];

    % Geometry
    opt.columnLength        = 53.6e-2;        % m
    opt.columnDiameter      = 2.60e-2;        % m
    opt.particleRadius      = 0.325e-4 /2;    % m
    opt.porosityColumn      = 0.38;
    opt.porosityParticle    = 0.00001;        % e_p very small to ensure e_t = e_c

    % Parameter units transformation
    % The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;   % m^2
    flowRate.recycle    = 0.1395e-6;      % m^3/s
    flowRate.desorbent1 = 0.0414e-6;      % m^3/s
    flowRate.extract1   = 0.0348e-6;      % m^3/s
    flowRate.feed1      = 0.02e-6  ;      % m^3/s
    flowRate.raffinate1 = 0.0266e-6;      % m^3/s
    flowRate.desorbent2 = 0.0414e-6;      % m^3/s
    flowRate.extract2   = 0.0348e-6;      % m^3/s
    flowRate.feed2      = flowRate.raffinate1; % m^3/s
    flowRate.raffinate2 = 0.0332e-6;      % m^3/s
    opt.flowRate_extract1   = flowRate.extract1;
    opt.flowRate_extract2   = flowRate.extract2;
    opt.flowRate_raffinate2 = flowRate.raffinate2;

    % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);       % m/s
    interstVelocity.desorbent1= flowRate.desorbent1 / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract1  = flowRate.extract1 / (crossArea*opt.porosityColumn);      % m/s
    interstVelocity.feed1     = flowRate.feed1 / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate1= flowRate.raffinate1 / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent2= flowRate.desorbent2 / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract2  = flowRate.extract2 / (crossArea*opt.porosityColumn);      % m/s
    interstVelocity.feed2     = flowRate.feed2 / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate2= flowRate.raffinate2 / (crossArea*opt.porosityColumn);    % m/s

    concentrationFeed 	= [0.5, 0.5, 0.5];   % g/m^3 [concentration_compA, concentration_compB]
    opt.molMass         = [180.16, 180.16, 180.16];
    opt.yLim            = max(concentrationFeed ./ opt.molMass);

    % Feed concentration setup
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), opt.nComponents);
    Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);
    for i = 1:opt.nComponents
        Feed.concentration(1:end, i) = concentrationFeed(i) / opt.molMass(i);
    end

    if opt.nZone == 8
        Desorbent = cell(1,2);
        Desorbent{1}.concentration = zeros(length(Feed.time), opt.nComponents);
        Desorbent{2}.concentration = zeros(length(Feed.time), opt.nComponents);
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
