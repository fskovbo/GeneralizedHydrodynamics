classdef LiebLinigerSolver < GeneralizedHydroSolver
    % Class specifying model specific quantities.
    % Solves the Bethe-Boltzmann equation for GHD propagation via the
    % superclass, which implements the general methods. 
    %
    % The units used in this parameterization are as follows:
    % m         = 1/2
    % hbar      = 1
    % g_1d      = 1
    % Lg        = hbar^2/(m*g_1d) = 1 (unit of length)
    % Eg        = 0.5*m*g_1d^2/hbar^2 = 1 (unit of energy)
    % rapidity  = k, whereby p = hbar*k = rapid
    %  
    
properties (Access = protected)
    
    % Formulas of model formated as strings
    energy      = 'rapid^2 - mu';
    momentum    = 'rapid';
    scattering  = '-2*atan(rapid/c)'
    
    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
    % Names of couplings in cell (must match model formulas)
    couplingNames = {'mu' , 'c'};
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = LiebLinigerSolver(x_grid, rapid_grid, couplings, stepOrder, extrapFlag)   
        % Set default parameters
        if nargin < 4
            % Propagate theta using stepper of order X.
            stepOrder = 2;
        end
        if nargin < 5
            % Extrapolate theta when propagating (necessary for homogeneous states)
            extrapFlag = false;
        end
        
        % Lieb-Liniger model has 1 species of quasi-particles
        Ntypes = 1;
   
        % Call superclass constructor
        obj = obj@GeneralizedHydroSolver(x_grid, rapid_grid, couplings, Ntypes, stepOrder, extrapFlag);
    end
    
    
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, setCouplingFlag)
        % Finds mu_0 for a given potential, V_ext, and temperature, T,
        % corresponding to a given atomnumber.
        % NOTE: V_ext is anonymous function with argument (t,x).
        
        if nargin < 5
            setCouplingFlag = false;
        end
        
        % Fit mu0 to Natoms
        mu0_guess   = 0;
        fitfunc     = @(mu0) abs( Natoms - calcNA(obj, mu0, T, V_ext) );
        mu0_fit     = fminsearch(fitfunc, mu0_guess);
        
        if setCouplingFlag % adjust couplings to result
            couplings_new   = obj.getCouplings();
            couplings_new.mu= @(t,x) mu0_fit - V_ext(t,x); 
            obj.setCouplings(couplings_new);
        end
        
        function Natoms_fit = calcNA(obj, mu0, T, V_ext)
            % Calculates number of atoms in stationary TBA state given by
            % specied paramters.
            couplings_fit   = obj.getCouplings();
            couplings_fit.mu= @(t,x) mu0 - V_ext(t,x);
            theta           = obj.calcThermalState(T, couplings_fit);
            density         = obj.calcCharges(theta, 0, 0);
            Natoms_fit      = trapz(permute(obj.x_grid, [5 2 3 4 1]), density);
        end % end nested function
    end
    
      
end % end public methods


methods (Access = protected)

    function h_i = getOneParticleEV(obj, t, x, rapid, i)
        h_i = rapid.^i;
    end
    

end % end private methods
    
end % end classdef