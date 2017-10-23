function [opt, interstVelocity, Feed, Desorbent] = getParameters(varargin)
% =============================================================================
%   Case 4, a 4-zone eight-column case for binary separation with SMA isotherm
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
    opt.nCellsParticle  = 5;
    opt.ABSTOL          = 1e-10;
    opt.INIT_STEP_SIZE  = 1e-14;
    opt.MAX_STEPS       = 5e6;

    % The parameter setting for the SMB
    opt.nInterval       = 5;
    opt.switch          = 257/opt.nInterval; % s
    opt.timePoints      = 1000/opt.nInterval;
    opt.Purity_limit    = [0.99, 0.99];
    opt.Penalty_factor  = 10;

    opt.enableDebug = true;
    opt.nZone       = 4;
    opt.nColumn     = 8;
    opt.structID    = [2 2 2 2]; % the column configuration which is used for structure optimization

    % Binding: Linear Binding isotherm
    opt.BindingModel = 'StericMassActionBinding';
    opt.nComponents = 2;
    opt.LAMBDA = 210; % Ionic capacity in mol/m^3
    opt.KA = [0, 0.1792, 10.83]; % [MB, BSA], MB for raffinate, BSA for extract, component 1 is salt
    opt.KD = [0, 1, 1]; % component 1 is salt
    opt.NU = [1, 1.22, 6.03]; % component 1 is salt
    opt.SIGMA = [1, 10, 75]; % component 1 is salt
    opt.compTargID = [2, 1]; % target compoents at [Extract Raffinate] ports

    % Transport
    opt.dispersionColumn          = [2.808, 1.755, 3.685, 1.755] .* 1e-7; % D_{ax}
    opt.filmDiffusion             = [1e-6, 1e-6, 1e-6];      % K_f
    opt.diffusionParticle         = [1.5, 3.6, 1.57] .* 1e-11;  % D_p
    opt.diffusionParticleSurface  = [0.0, 0.0, 0.0];

    % Geometry
    opt.columnLength        = 9.04e-2; % m
    opt.columnDiameter      = 1.0e-2;  % m
    opt.particleRadius      = 5.0e-5;  % m 
    opt.porosityColumn      = 0.39;
    opt.porosityParticle    = 0.57;

    % Parameter units transformation
    % The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;   % m^2
    flowRate.recycle    = 5.4812e-8;      % m^3/s
    flowRate.feed       = 1.03e-6/60;      % m^3/s
    flowRate.raffinate  = 1.06e-6/60;      % m^3/s
    flowRate.desorbent  = 2.05e-6/60;      % m^3/s
    flowRate.extract    = 2.02e-6/60;      % m^3/s
    opt.flowRate_extract   = flowRate.extract;
    opt.flowRate_raffinate = flowRate.raffinate;

    % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s

    concentrationFeed 	= [0.1, 0.5];   % g/m^3 [MB BSA]
    opt.molMass         = [66.463, 66.463]; % The molar mass of each components
    opt.yLim            = max(concentrationFeed ./ opt.molMass) * 1.1; % the magnitude for plotting

    % Feed concentration setup
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), opt.nComponents);
    Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);

    for i = 1:opt.nComponents
        Feed.concentration(1:end, i) = concentrationFeed(i) / opt.molMass(i);
    end

    % With steric mass action isotherm, salt concentration should be firstly appended
    if strcmp('StericMassActionBinding', opt.BindingModel)
        opt.concentrationSalt  = [130, 340]; % mol/m^3
        Feed.concentration = [ones(length(Feed.time), 1).*opt.concentrationSalt(1),  Feed.concentration];
        Desorbent.concentration = [ones(length(Feed.time), 1).*opt.concentrationSalt(2), Desorbent.concentration];
    end

% -----------------------------------------------------------------------------
    % Capable of placing a CSTR or DPFR apparatues before and after the calculated column

    % Continuous Stirred Tank Reactor
    opt.enable_CSTR = false;
    opt.CSTR_length = 0.01;

    % Dispersive Plug Flow Reactor
    opt.enable_DPFR = false;

    opt.DPFR_length = 0.0016;
    opt.DPFR_nCells = 50;

    opt.DPFR_velocity   = 0.00315;
    opt.DPFR_dispersion = 2.5e-30;


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
