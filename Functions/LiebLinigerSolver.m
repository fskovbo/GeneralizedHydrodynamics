classdef LiebLinigerSolver < GeneralizedHydroSolver
    % Class specifying model specific quantities.
    % Solves the Bethe-Boltzmann equation for GHD propagation via the
    % superclass, which implements the general methods. 
    
properties (Access = private)
    
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
   
        % Call superclass constructor
        obj = obj@GeneralizedHydroSolver(x_grid, rapid_grid, couplings, stepOrder, extrapFlag);
    end
    
    
    function setCouplings(obj, couplings)
        % Check that couplings match Lieb-Liniger model
        % Lieb-Liniger couplings are chemical potential, mu, and particle
        % interaction, c, along with their derivatives.
        LLcouplings = {'mu', 'dmudt', 'dmudx', 'c', 'dcdt', 'dcdx'};
        fn = fieldnames(couplings);
        for i = 1:length(LLcouplings)                      
            if ~any(strcmp(fn, LLcouplings{i}))
                error([LLcouplings{i} ' is not defined in couplings struct!'])
            end
        end
        
        obj.couplings   = couplings;
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
            couplings_fit   = obj.getCouplings();
            couplings_fit.mu= @(t,x) mu0 - V_ext(t,x);
            theta           = obj.calcThermalState(T, couplings_fit);
            density         = obj.calcCharges(theta, 0, 0);
            Natoms_fit      = trapz(permute(obj.x_grid, [3 2 1]), density);
        end % end nested function
    end
    
      
end % end public methods


methods (Access = protected)
    
    function ebare = getBareEnergy(obj, t, x, rapid)
        % m = 1/2
        ebare = rapid.^2 - obj.couplings.mu(t,x);
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid)
        pbare = rapid;
    end
    
    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2)
        % for now: rapid1 is 2nd dim, rapid2 is 1st dim
        if size(rapid1,1) ~= 1 && numel(rapid1) ~= 1
            rapid1 = permute(rapid1, [2 1 3]);
        end
        if size(rapid2,2) ~= 1 && numel(rapid2) ~= 1
            rapid2 = permute(rapid2, [2 1 3]);
        end
        
        dT = -2*obj.couplings.c(t,x)./(obj.couplings.c(t,x).^2 + (rapid1 - rapid2).^2);
        dT(isnan(dT)) = 0; % removes any NaN
    end
    
    
    function de = calcEnergyRapidDeriv(obj, t, x, rapid)
        % m = 1/2
        de = 2*rapid;
    end
    
    
    function dp = calcMomentumRapidDeriv(obj, t, x, rapid)
        dp = ones(1,obj.N);
    end
    
    
    function f = getStatFactor(obj, theta)
        % Quasi-particles in LL-model are fermions
        f = 1 - theta;
    end
    
    
    function h_i = getOneParticleEV(obj, i, rapid)
        % NOTE: This might not be correct (doesnt seem right for i = 2 (energy))
        h_i = rapid.^i /factorial(i);
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
       theta = 1./(1 + exp(e_eff)); 
    end
    
    
    function [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid)        
        % Calculates velocities, outputs (1xNxM)
        de_dr   = obj.applyDressing(obj.calcEnergyRapidDeriv(t, x, rapid), theta, t);
        dp_dr   = obj.applyDressing(obj.calcMomentumRapidDeriv(t, x, rapid), theta, t);
        
        % Calculate acceleration from inhomogenous potential. Note dmudt
        % does not contribute as f = 0;
        a_eff_mu = obj.couplings.dmudx(t,x);
        
        if size(a_eff_mu,2) == 1
            a_eff_mu = repmat(a_eff_mu, 1, obj.N); % should be (1xNxM)
        end
        
        % Calculate acceleration from inhomogenous interaction
        % NEED TO CHECK DIMENSIONS/INDICES OF RAPID!!
        a_eff_c = 0;
        if ~isempty(obj.couplings.dcdt) || ~isempty(obj.couplings.dcdx)
            % Calculate derivative of scattering phase with respect to
            % interaction c
            dk      = obj.rapid_grid(2) - obj.rapid_grid(1);            
            dTdc    = -2*(rapid - obj.rapid_grid')./(obj.couplings.c(t,x).^2 + (rapid - obj.rapid_grid').^2);
            dTdc(isnan(dTdc)) = 0; % removes any NaN
            B       = dk/(2*pi) * dTdc.*theta;
        end
        
        if ~isempty(obj.couplings.dcdt) % calc time deriv contribution
            f = zeros(obj.N, 1, obj.M);
            for i = 1:obj.M
                % Transpose to column vector for matrix multiplication
                f(:,:,i) = B(:,:,i) * transpose(dp_dr(:,:,i));
            end
        
            f_dr    = obj.applyDressing(permute(f, [2 1 3]), theta, t);
            a_eff_c = a_eff_c + obj.couplings.dcdt(t,x).*f_dr;
        end

        if ~isempty(obj.couplings.dcdx) % calc space deriv contribution
            L = zeros(obj.N, 1, obj.M);
            for i = 1:obj.M
                % Transpose to column vector for matrix multiplication
                L(:,:,i) = B(:,:,i) * transpose(de_dr(:,:,i));
            end
        
            L_dr    = obj.applyDressing(permute(L, [2 1 3]), theta, t);
            a_eff_c = a_eff_c + obj.couplings.dcdx(t,x).*L_dr;
        end
        
        a_eff_c = a_eff_c./dp_dr;
        
        v_eff   = de_dr./dp_dr;
        a_eff   = a_eff_c + a_eff_mu;
    end
    
    
end % end private methods
    
end % end classdef