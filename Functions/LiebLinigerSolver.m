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
    function obj = LiebLinigerSolver(x_grid, rapid_grid, couplings, Options)   
        % Set default parameters
        if nargin < 4
            Options = struct;
        end
        
        % Lieb-Liniger model has 1 species of quasi-particles
        Ntypes = 1;
   
        % Call superclass constructor
        obj = obj@GeneralizedHydroSolver(x_grid, rapid_grid, couplings, Ntypes, Options);
    end
    
    
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, setCouplingFlag)
        % Finds mu_0 for a given potential, V_ext, and temperature, T,
        % corresponding to a given atomnumber.
        % NOTE: V_ext is anonymous function with argument (t,x).
        
        if nargin < 5
            setCouplingFlag = false;
        end
        
        if isempty(V_ext)
            V_ext = obj.couplings{1};
        end
        
        % Fit mu0 to Natoms
        mu0_guess   = 0;
        fitfunc     = @(mu0) abs( Natoms - calcNA(obj, mu0, T, V_ext) );
        options     = optimset('Display','iter');
        mu0_fit     = fminsearch(fitfunc, mu0_guess,options);
        
        if setCouplingFlag % adjust couplings to result
            couplings_new   = obj.getCouplings();
            couplings_new{1}= @(t,x) mu0_fit - V_ext(t,x);

            obj.setCouplings(couplings_new);
        end
        
        function Natoms_fit = calcNA(obj, mu0, T, V_ext)
            % Calculates number of atoms in stationary TBA state given by
            % specied paramters.
            couplings_fit   = obj.getCouplings();
            couplings_fit{1}= @(t,x) mu0 - V_ext(t,x);
            theta           = obj.calcThermalState(T, couplings_fit);
            density         = obj.calcCharges(theta, 0, 0);
            Natoms_fit      = trapz(permute(obj.x_grid, [5 2 3 4 1]), density);
        end % end nested function
    end
    
    
    function nk = calcMomentumDistr(obj, theta, t_array)
        
        Nsteps  = length(theta); % number of time steps
        nk       = zeros(obj.N, Nsteps);
        delta_x = obj.x_grid(2) - obj.x_grid(1);
        rho     = obj.transform2rho(theta, t_array);
        
        if ~iscell(rho)
            nk(:,1) = delta_x * squeeze(sum( double(rho) , 5 ));
        else
            for i = 1:Nsteps
                rho_i = rho{i};
                nk(:,i) = delta_x * squeeze(sum( double(rho_i) , 5 ));
            end
        end
    end
    
      
end % end public methods


methods (Access = protected)

    function h_i = getOneParticleEV(obj, t, x, rapid, i)
        h_i = rapid.^i;
    end
    
    
    function [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid, type)        
        % Overloads method, as acceleration from mu has much simpler
        % expression than the general one.
        de_dr   = obj.applyDressing(obj.calcEnergyRapidDeriv(t, x, rapid, type), theta, t);
        dp_dr   = obj.applyDressing(obj.calcMomentumRapidDeriv(t, x, rapid, type), theta, t);
        
        v_eff   = de_dr./dp_dr;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = GHDtensor( zeros(size( v_eff )) );
            return
        end
        
        % Calculate acceleration from inhomogenous potential. Note dmudt
        % does not contribute as f = 0;
        a_eff_mu = 0;
        if ~isempty(obj.couplingDerivs{2,1})
            a_eff_mu = obj.couplingDerivs{2,1}(t,x);
            if size(a_eff_mu,1) == 1
                a_eff_mu = repmat(a_eff_mu, obj.N, 1); 
            end
        end
        a_eff_mu = GHDtensor(a_eff_mu);
        
        % Calculate acceleration from inhomogenous interaction
        a_eff_c = 0;
        if ~isempty(obj.couplingDerivs{1,2}) || ~isempty(obj.couplingDerivs{2,2})
            % Calculate derivative of scattering phase with respect to
            % interaction c
            delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);            
            dT      = obj.calcScatteringCouplingDeriv(2, t, x, rapid, obj.rapid_grid, type, obj.type_grid);
            B       = delta_k/(2*pi) * dT.*transpose(theta);
        end
        
        if ~isempty(obj.couplingDerivs{1,2}) % calc time deriv contribution
            f       = B*dp_dr;
            f_dr    = obj.applyDressing(f, theta, t);
            a_eff_c = a_eff_c + obj.couplingDerivs{1,2}(t,x).*f_dr;
        end

        if ~isempty(obj.couplingDerivs{2,2}) % calc space deriv contribution
            L       = B*de_dr;
            L_dr    = obj.applyDressing(L, theta, t);
            a_eff_c = a_eff_c + obj.couplingDerivs{2,2}(t,x).*L_dr;
        end
        
        a_eff_c = a_eff_c./dp_dr;
        a_eff   = a_eff_c + a_eff_mu;
    end
    

end % end private methods
    
end % end classdef