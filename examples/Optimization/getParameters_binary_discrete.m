function [opt, interstVelocity, Feed] = getParameters(ParSwarm)
%   Case 1, a four-column demonstration case

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


    valueAssign = struct('columnLength',ParSwarm(1), 'switch',ParSwarm(2),...
        'recycle',ParSwarm(3),'feed',ParSwarm(4), 'desorbent',ParSwarm(5),...
        'extract',ParSwarm(6));

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
    opt.switch          = valueAssign.switch;  % switching time 
    opt.timePoints      = 1000;         % the observed time-points
    opt.Purity_extract_limit   = 0.9964;  % used for constructing constraints
    opt.Purity_raffinate_limit = 0.9975;  % used for constructing constraints
    opt.Penalty_factor         = 10;    % penalty factor in penalty function

    opt.enableDebug = false; % set it false if you are using the optimizer
    opt.nZone   = 4;
    opt.nColumn = 6; % The amount of columns are known
    opt.structNumber = 20;
%     opt.structID = [1 1 1 1]; % the opt.structID is unknown

%   Binding: Linear Binding isotherm
    opt.BindingModel = 'LinearBinding';
    opt.nComponents = 2;
    opt.KA = [5.72 7.7]; % [comp_A, comp_B], A for raffinate, B for extract
    opt.KD = [1, 1];
    opt.comp_raf_ID = 1; % the target component withdrawn from the raffinate ports
    opt.comp_ext_ID = 2; % the target component withdrawn from the extract ports

%   Transport
    opt.dispersionColumn          = 3.8148e-20;     % D_{ax}
    opt.filmDiffusion             = [100 100];      % K_{eff} 
    opt.diffusionParticle         = [1.6e4 1.6e4];  % D_p
    opt.diffusionParticleSurface  = [0.0 0.0];

%   Geometry
    opt.columnLength        = valueAssign.columnLength;      % m
    opt.columnDiameter      = 0.02;      % m
    opt.particleRadius      = 0.0005;    % m % user-defined one in this case
    opt.porosityColumn      = 0.83;
    opt.porosityParticle    = 0.000001;  % unknown

%   Parameter units transformation
%   The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;      % m^2
    flowRate.recycle    = valueAssign.recycle;      % m^3/s  
    flowRate.feed       = valueAssign.feed;         % m^3/s
    flowRate.desorbent  = valueAssign.desorbent;    % m^3/s
    flowRate.extract    = valueAssign.extract;      % m^3/s
    flowRate.raffinate  = flowRate.desorbent - flowRate.extract + flowRate.feed;        % m^3/s
    opt.flowRate_extract   = flowRate.extract;
    opt.flowRate_raffinate = flowRate.raffinate;

%   Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s 
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s

    concentrationFeed 	= [0.55, 0.55];   % g/m^3 [concentration_compA, concentration_compB]
    opt.molMass         = [180.16, 180.16]; % The molar mass of each components
    opt.yLim            = max(concentrationFeed ./ opt.molMass); % the magnitude for plotting

%   Feed concentration setup   
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), opt.nComponents);

    for i = 1:opt.nComponents
        Feed.concentration(1:end,i) = (concentrationFeed(i) / opt.molMass(i));
    end

% -----------------------------------------------------------------------------
%   Capable of placing a CSTR or DPFR apparatues before and after the calculated column

%   Continuous Stirred Tank Reactor
    opt.enable_CSTR = false;
    opt.CSTR_length = 0.01;

%   Dispersive Plug Flow Reactor
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
%      Copyright © 2008-2016: Eric von Lieres, Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
