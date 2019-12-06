classdef sinhGordonSolver < GeneralizedHydroSolver
    % Class specifying model specific quantities.
    % Solves the Bethe-Boltzmann equation for GHD propagation via the
    % superclass, which implements the general methods. 

    
properties (Access = protected)
    
    % Formulas of model formated as strings
    energy      = 'cosh(rapid) - mu';
    momentum    = 'sinh(rapid)';
    scattering  = '1i*log( tanh(1/2 * (rapid - 1i*pi*B/2)) / tanh(1/2 * (rapid + 1i*pi*B/2)) )'
    
    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
    % Names of couplings in cell (must match model formulas)
    couplingNames = {'mu' , 'B'};
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = sinhGordonSolver(x_grid, rapid_grid, rapid_w, couplings, Options)   
        % Set default parameters
        if nargin < 5
            Options = struct;
        end
        
        % Lieb-Liniger model has 1 species of quasi-particles
        Ntypes = 1;
   
        % Call superclass constructor
        obj = obj@GeneralizedHydroSolver(x_grid, rapid_grid, rapid_w, couplings, Ntypes, Options);
    end
    
      
end % end public methods


methods (Access = protected)

    function h_i = getOneParticleEV(obj, t, x, rapid, i)
        switch i
        case 0 % eigenvalue of S_z operator
            h_i = repmat(obj.type_grid, obj.N, 1);
        case 1 
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of H operator WITHOUT B-field
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid);
        otherwise 
            error(['Eigenvalue ' num2str(i) ' not implmented!'])
        end
    end
    

end % end private methods
    
end % end classdef