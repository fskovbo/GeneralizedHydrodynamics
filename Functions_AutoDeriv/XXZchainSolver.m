classdef XXZchainSolver < GeneralizedHydroSolver
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
    energy      = '-sinh(Theta) * sinh(type*Theta) / (cosh(type*Theta) - cos(2*rapid)) - type*B';
    momentum    = '2*atan( coth(type*Theta/2) * tan(rapid) )';
    scattering  = '0' % overloaded derivatives
    
    % Species of quasiparticle
    quasiSpecies= 'fermion'; 
    
    % Names of couplings in cell (must match model formulas)
    couplingNames = {'B' , 'Theta'};
    
end % end private properties
    
    
methods (Access = public)
        
    % Constructor
    function obj = XXZchainSolver(x_grid, rapid_grid, couplings, Ntypes, stepOrder, extrapFlag)   
        % Set default parameters
        if nargin < 5
            % Propagate theta using stepper of order X.
            stepOrder = 2;
        end
        if nargin < 6
            % Extrapolate theta when propagating (necessary for homogeneous states)
            extrapFlag = false;
        end
   
        % Call superclass constructor
        obj = obj@GeneralizedHydroSolver(x_grid, rapid_grid, couplings, Ntypes, stepOrder, extrapFlag);
        
        % XXZ model has in principle infinite species of quasi-particles
    end
    
      
end % end public methods


methods (Access = protected)

    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.calcMomentumRapidDeriv(t, x, r_arg, abs(type1-type2));
        dT2     = obj.calcMomentumRapidDeriv(t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.calcMomentumRapidDeriv(t, x, r_arg, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = GHDtensor( dT );
    end
    

    function dT = calcScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        if coupIdx == 1 % deriv w.r.t. B
            dT  = 0;
            dT  = GHDtensor( dT );
            return
        end 
        % Else it's deriv w.r.t. Delta
        
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.calcMomentumCouplingDeriv(2, t, x, r_arg, abs(type1-type2));
        dT2     = obj.calcMomentumCouplingDeriv(2, t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.calcMomentumCouplingDeriv(2, t, x, r_arg, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = GHDtensor( dT );
    end
    
    
    function h_i = getOneParticleEV(obj, t, x, rapid, i)
        switch i
        case 0 % eigenvalue of S_z operator
            h_i = repmat(obj.type_grid, obj.N, 1);
        case 1 
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of H operator
            h_i = obj.getBareEnergy(t, x, rapid,  obj.type_grid) + obj.type_grid.*obj.couplings{1}(t,x);
        otherwise 
            error(['Eigenvalue ' num2str(i) ' not implmented!'])
        end
    end
    

    
    
end % end private methods
    
end % end classdef