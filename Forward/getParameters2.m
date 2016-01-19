function [opt, interstVelocity, Feed] = getParameters()

% =============================================================================
% This is the function to input all the necessary data for simulation

% Returns: 
%       1. opt stands for options, which involves the parameter settings
%       for the algorithm, the binding isotherm, and the model equations

%       2. interstVelocity is calculated from flowrate of each column and inlet. 
%       interstitial_velocity = flow_rate / (across_area * porosity_Column)

%       3. Feed initializes the injection concentration
% =============================================================================


%   the parameter setting for simulator
    opt.tolIter         = 1e-4;
    opt.nMaxIter        = 1000;
    opt.nThreads        = 8;
    opt.nCellsColumn    = 40;
    opt.nCellsParticle  = 1;
    opt.switch          = 1552;
    opt.timePoints      = 1000;
    opt.ABSTOL          = 1e-10;
    opt.INIT_STEP_SIZE  = 1e-14;
    opt.MAX_STEPS       = 5e6;
    
    opt.Purity_extract_limit    = 0.99;
    opt.Purity_raffinate_limit  = 0.99;
    opt.Penalty_factor          = 10;

    opt.enableDebug = true;
    opt.nColumn = 4; 
%     opt.nColumn = 8;

%   Binding: Linear Binding isotherm
    opt.nComponents = 2;
    opt.KA = [0.28 0.54];
    opt.KD = [1, 1];
    
%   Transport
    opt.dispersionColumn          = 3.8148e-6;      % 
    opt.filmDiffusion             = [100 100];      % unknown 
    opt.diffusionParticle         = [1.6e4 1.6e4];  % unknown
    opt.diffusionParticleSurface  = [0.0 0.0];

%   Geometry
    opt.columnLength        = 53.6e-2;        % m
    opt.columnDiameter      = 2.60e-2;        % m
    opt.particleRadius      = 0.325e-2 /2;    % m
    opt.porosityColumn      = 0.38;
    opt.porosityParticle    = 0.00001;        % unknown

%   Parameter units transformation
%   The flow rate of Zone I was defined as the recycle flow rate
    crossArea = pi * (opt.columnDiameter/2)^2;        % m^2
    flowRate.recycle    = 0.1395e-6;      % m^3/s 
    flowRate.feed       = 0.02e-6  ;      % m^3/s
    flowRate.raffinate  = 0.0266e-6;      % m^3/s
    flowRate.desorbent  = 0.0414e-6;      % m^3/s
    flowRate.extract    = 0.0348e-6;      % m^3/s
    
%   Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
    interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s 
    interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
    interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
    interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
   
    concentrationFeed 	= [0.5, 0.5];   % g/m^3 [concentration_compA, concentration_compB]
    opt.FructoseMolMass = 262.1535; % g/mol
    opt.GlucoseMolMass  = 262.1535; % g/mol
    opt.flowRate_recycle = flowRate.recycle;
    
    
%   Feed concentration setup    
    Feed.time = linspace(0, opt.switch, opt.timePoints);
    Feed.concentration = zeros(length(Feed.time), 2);

    Feed.concentration(1:end,1) = (concentrationFeed(1) / opt.FructoseMolMass );
    Feed.concentration(1:end,2) = (concentrationFeed(2) / opt.GlucoseMolMass);

end
% =============================================================================
%  SMB - The Simulated Moving Bed Chromatography for separation of
%  target compounds, such as fructose and glucose.
%  
%  Author: QiaoLe He   E-mail: q.he@fz-juelich.de
%                                      
%  Institute: Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%  
%  All rights reserved. Please see the license of CADET.
% =============================================================================
