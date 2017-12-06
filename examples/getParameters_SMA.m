function [opt, interstVelocity, Feed, Desorbent] = getParameters(iter, varargin)
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
    opt.switch          = 160; % s
    opt.timePoints      = 1000;
    opt.Purity_limit    = [0.99, 0.99];
    opt.Penalty_factor  = 10;

    opt.enableDebug = true;
    opt.nZone       = 4;
    opt.nColumn     = 8;
    opt.structID    = ones(1, opt.nZone) .* 2;

    % Binding: Linear Binding isotherm
    opt.BindingModel = 'StericMassActionBinding';

    % Geometry
    opt.columnLength        = 2.4e-2; % m
    opt.columnDiameter      = 1e-2;   % m
    opt.particleRadius      = 4.5e-5; % m
    opt.porosityColumn      = 0.37;
    opt.porosityParticle    = 0.75;

    if iter == 1

        opt.nComponents = 3;
        opt.LAMBDA = 1.2e3; % mol/m^3
        opt.KA = [0, 1.59, 7.70, 35.5];
        opt.KD = [0, 1000, 1000, 1000];
        opt.NU = [1, 5.29, 3.70, 4.70];
        opt.SIGMA = [1, 10.60, 10.0, 11.83];
        opt.compTargID = [3, 1];

        % Transport
        opt.dispersionColumn          = ones(1, opt.nZone) .* 5.75e-8; % D_{ax}
        opt.filmDiffusion             = ones(1, opt.nComponents+1) .* 6.90e-6;   % K_f
        opt.diffusionParticle         = ones(1, opt.nComponents+1) .* 6.07e-11;  % D_p
        opt.diffusionParticleSurface  = zeros(1, opt.nComponents+1);

        % Parameter units transformation
        % The flow rate of Zone I was defined as the recycle flow rate
        crossArea = pi * (opt.columnDiameter/2)^2;   % m^2
        flowRate.recycle   = 3.21e-8;      % m^3/s
        flowRate.desorbent = 2.39e-8;      % m^3/s
        flowRate.extract   = 1.09e-8;      % m^3/s
        flowRate.feed      = 0.55e-8;      % m^3/s
        flowRate.raffinate = 0.85e-8;      % m^3/s
        opt.flowRate_raffinate = flowRate.raffinate;
        opt.flowRate_extract   = flowRate.extract;

        % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
        interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
        interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
        interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
        interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
        interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s

        concentrationFeed 	= [0.1, 0.1, 0.1];   % g/m^3 [MB BSA]
        opt.molMass         = [66.463, 66.463, 66.463];
        opt.yLim            = max(concentrationFeed ./ opt.molMass) .* 1.1;

        % Feed and Desorbent concentration setup
        Feed.time = linspace(0, opt.switch, opt.timePoints);
        Feed.concentration = zeros(length(Feed.time), opt.nComponents);
        Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);

        for i = 1:opt.nComponents
            Feed.concentration(1:end, i) = concentrationFeed(i) / opt.molMass(i);
        end

        % With steric mass action isotherm, salt concentration should be firstly appended
        if strcmp('StericMassActionBinding', opt.BindingModel)
            opt.concentrationSalt  = [100, 370]; % mol/m^3
            Feed.concentration = [ones(length(Feed.time), 1).*opt.concentrationSalt(1),  Feed.concentration];
            Desorbent.concentration = [ones(length(Feed.time), 1).*opt.concentrationSalt(2), Desorbent.concentration];
        end

    elseif iter == 2

        opt.nComponents = 2;
        opt.LAMBDA = 1.2e3; % mol/m^3
        opt.KA = [0, 1.59, 7.70];
        opt.KD = [0, 1000, 1000];
        opt.NU = [1, 5.29, 3.70];
        opt.SIGMA = [1, 10.60, 10.0];
        opt.compTargID = [1, 2];

        % Transport
        opt.dispersionColumn          = ones(1, opt.nZone) .* 5.75e-8; % D_{ax}
        opt.filmDiffusion             = ones(1, opt.nComponents+1) .* 6.90e-6;   % K_f
        opt.diffusionParticle         = ones(1, opt.nComponents+1) .* 6.07e-11;  % D_p
        opt.diffusionParticleSurface  = zeros(1, opt.nComponents+1);

        % Parameter units transformation
        % The flow rate of Zone I was defined as the recycle flow rate
        crossArea = pi * (opt.columnDiameter/2)^2;   % m^2
        flowRate.recycle   = 2.56e-8;      % m^3/s
        flowRate.desorbent = 1.24e-8;      % m^3/s
        flowRate.extract   = 1.09e-8;      % m^3/s
        flowRate.feed      = 0.60e-8;      % m^3/s % keep it same as the flowRate.raffinate when iter = 1
        flowRate.raffinate = 0.75e-8;      % m^3/s
        opt.flowRate_raffinate = flowRate.raffinate;
        opt.flowRate_extract   = flowRate.extract;

        % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
        interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
        interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
        interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
        interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
        interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s

        concentrationFeed 	= [0.1, 0.1];   % g/m^3 [MB BSA]
        opt.molMass         = [66.463, 66.463];
        opt.yLim            = max(concentrationFeed ./ opt.molMass) .* 1.5;

        % Feed and Desorbent concentration setup
        Feed.time = linspace(0, opt.switch, opt.timePoints);
        Feed.concentration = zeros(length(Feed.time), opt.nComponents);
        Desorbent.concentration = zeros(length(Feed.time),  opt.nComponents);

        % With steric mass action isotherm, salt concentration should be firstly appended
        if strcmp('StericMassActionBinding', opt.BindingModel)
            opt.concentrationSalt  = [200, 260]; % mol/m^3
            Desorbent.concentration = [ones(length(Feed.time), 1).* opt.concentrationSalt(2), Desorbent.concentration];
        end

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
