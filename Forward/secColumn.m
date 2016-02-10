function [outletProfile, lastState] = secColumn(inletProfile, params, lastState)

% =============================================================================
% Simulation of the single column
%
% Parameters:
%       - inletProfile. Inlet time and corresponding concentration
%       - params. Get parameters for simulation
%       - lastState. The recorded last STATE from previous simulation
%       of next simulation
% 
% Returns:
%       - outletProfile. outlet time and corresponding concentration
%       - lastState. Record the last STATE which used as the initial state
% =============================================================================


    if nargin < 3
        lastState = [];
    end
    
    if isempty(params.initMobilCon) && isempty(params.initSolidCon) && isempty(lastState)
        warning('There are no Initial Conditions / Boundary Conditions for the Simulator');
    end
    
%   Get parameters
    [opt, ~, ~] = getParameters();
    
    model = ModelGRM();
    model.nComponents = opt.nComponents;
    model.kineticBindingModel = false;
    model.bindingModel = LinearBinding();
       
%   Adsorption parameters
    model.bindingParameters.LIN_KA         = opt.KA;
    model.bindingParameters.LIN_KD         = opt.KD;
    
    if nargin >= 3 && ~isempty(lastState)
        model.initialState = lastState;
    else      
        model.initialMobileConcentration = params.initMobilCon;
        model.initialSolidConcentration  = params.initSolidCon;
    end
    
%   Transport
    model.dispersionColumn          = opt.dispersionColumn;
    model.filmDiffusion             = opt.filmDiffusion;
    model.diffusionParticle         = opt.diffusionParticle;
    model.diffusionParticleSurface  = opt.diffusionParticleSurface;
    model.interstitialVelocity      = params.interstitialVelocity;

%   Geometry
    model.columnLength        = opt.columnLength;
    model.particleRadius      = opt.particleRadius;
    model.porosityColumn      = opt.porosityColumn;
    model.porosityParticle    = opt.porosityParticle;
    
%   Apply the inlet profile to the CADET model
    Time = repmat({inletProfile.time}, 1, opt.nComponents);
    if opt.nComponents == 2
        Profile = [{inletProfile.concentration(:,1)}, {inletProfile.concentration(:,2)}];
    elseif opt.nComponents == 3
        Profile = [{inletProfile.concentration(:,1)}, {inletProfile.concentration(:,2)}, {inletProfile.concentration(:,3)}];
    end
    
    model.setInletsFromData(Time, Profile);

%   Turn off the warnings of the interpolation
    warning('off', 'MATLAB:interp1:ppGriddedInterpolant');
	warning('off', 'MATLAB:interp1:UsePCHIP');

%   Discretization
    disc = DiscretizationGRM();
    disc.nCellsColumn   = opt.nCellsColumn;
    disc.nCellsParticle = opt.nCellsParticle;
    
%   Solving options
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
   
    
%   Run the simulation
    try
        result = sim.simulate();
    catch e
        % Something went wrong
        error('CADET:simulationFailed', 'Check your settings and try again.\n%s', e.message);
    end
    
%   Extract the outlet profile
    outletProfile.time = result.solution.time;
    outletProfile.concentration = result.solution.outlet(:,:);
    lastState =  result.solution.lastState;
 
    
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