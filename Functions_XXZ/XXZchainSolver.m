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
    
properties (Access = private)
    
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
        % and has periodic boundary conditions for rapidity
        obj.periodRapid = true;
        
    end
    
    
    function setCouplings(obj, couplings)
        % Check that couplings match Lieb-Liniger model
        % Lieb-Liniger couplings are chemical potential, mu, and particle
        % interaction, c, along with their derivatives.
        XXZcouplings = {'B', 'dBdt', 'dBdx', 'Delta', 'dDeltadt', 'dDeltadx'};
        fn = fieldnames(couplings);
        for i = 1:length(XXZcouplings)                      
            if ~any(strcmp(fn, XXZcouplings{i}))
                error([XXZcouplings{i} ' is not defined in couplings struct!'])
            end
        end
        
        obj.couplings   = couplings;
    end
    
      
end % end public methods


methods (Access = protected)
    
    function ebare = getBareEnergy(obj, t, x, rapid, type)
        phi     = acosh(obj.couplings.Delta(t,x));
        ebare   = -0.5*sinh(phi).*obj.calcMomentumRapidDeriv(t, x, rapid, type) - type.*obj.couplings.B(t,x);
    end
    
    
    function ebare = getBareEnergyNoField(obj, t, x, rapid, type)
        phi     = acosh(obj.couplings.Delta(t,x));
        ebare   = -0.5*sinh(phi).*obj.calcMomentumRapidDeriv(t, x, rapid, type);
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        phi     = acosh(obj.couplings.Delta(t,x));
        pbare   = 2*atan( coth(type.*phi/2) .* tan(rapid) );
    end
    
    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Reshape input to right dimensions
%         rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
%         rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
%         type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
%         type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        if size(rapid1,2) ~= 1  % rapid1 is 1st index
            rapid1 = permute(rapid1, [2 1 3 4 5]); 
        end
        if size(rapid2,1) ~= 1  % rapid2 is 2nd index
            rapid2 = permute(rapid2, [2 1 3 4 5]); 
        end
        if size(type1,4) ~= 1  % type1 is 3rd index
            type1 = permute(type1, [1 2 4 3 5]); 
        end
        if size(type2,3) ~= 1  % type2 is 4th index
            type2 = permute(type2, [1 2 4 3 5]); 
        end
        
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
            ri_idx = min( i, size(r_arg, 3) );
            
            for j = 1:obj.Ntypes
                rj_idx = min( j, size(r_arg, 4) );
                
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.calcMomentumRapidDeriv(t, x, r_arg(:,:,ri_idx,rj_idx), n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = GHDtensor( dT );
    end
    
    
    function de = calcEnergyRapidDeriv(obj, t, x, rapid, type)
        phi = acosh(obj.couplings.Delta(t,x));
        de  = 2*sinh(phi).*sin(2*rapid).*sinh(type.*phi)./( cosh(type.*phi) - cos(2*rapid) ).^2;
    end
    
    
    function dp = calcMomentumRapidDeriv(obj, t, x, rapid, type)
        phi = acosh(obj.couplings.Delta(t,x));
        dp  = 2*sinh(type.*phi)./( cosh(type.*phi) - cos(2*rapid) );
    end
    
    
    function f = getStatFactor(obj, theta)
        % Quasi-particles in XXZ-model are fermions
        f = 1 - theta;
    end
    
    
    function dp = calcMomentumDeltaDeriv(obj, t, x, rapid, type)
        phi = acosh(obj.couplings.Delta(t,x));
%         dp  = csch(phi).*type.*sin(2*rapid)./(cos(2*rapid)-cosh(type.*phi));
        dp = -(type.*tan(rapid).*(coth((type.*phi)/2).^2 - 1))./(sinh(phi).*(coth((type.*phi)/2).^2.*tan(rapid).^2 + 1));
    end
    
    
    function de = calcEnergyDeltaDeriv(obj, t, x, rapid, type)
        phi = acosh(obj.couplings.Delta(t,x));
%         X   = cos(2*rapid)-cosh(type.*phi);
%         de  = csch(phi).*( ( -cosh(rapid).*sinh(type.*phi) - type.*sinh(rapid).*cosh(type.*phi))./X + type.*sinh(rapid).*sinh(type.*phi).^2 ./X.^2 );
        de = -(type.*sinh(phi) - cos(2*rapid).*sinh(type.*phi).*cosh(phi) + cosh(type.*phi).*sinh(type.*phi).*cosh(phi) - type.*cos(2*rapid).*cosh(type.*phi).*sinh(phi))./(sinh(phi).*(cos(2*rapid) - cosh(type.*phi)).^2);
    end
    
    
    function dT = calcScatteringDeltaDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % Reshape input to right dimensions        
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        % Calculate contribution from 2 first terms
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        
        r_arg   = (rapid1-rapid2);

        dT1     = (1 - I_type).*obj.calcMomentumDeltaDeriv(t, x, r_arg, abs(type1-type2));
        dT2     = obj.calcMomentumDeltaDeriv(t, x, r_arg, type1+type2);
        
        dT1(isnan(dT1)) = 0; % removes any NaN
        dT2(isnan(dT2)) = 0; % removes any NaN
        
        
        dT3 = zeros(size(dT1));
        for i = 1:obj.Ntypes
            for j = 1:obj.Ntypes
                for n = (abs(i-j)+2):2:(i+j-2)
                    temp = 2*obj.calcMomentumDeltaDeriv(t, x, r_arg, n);
                    temp(isnan(temp)) = 0; % removes any NaN
                    dT3(:,:,i,j,:) = dT3(:,:,i,j,:) + temp;
                end
            end
        end
        
        
        dT = dT1 + dT2 + dT3;
        
        dT  = GHDtensor( dT );
    end
    
    
    function F = getFreeEnergy(obj, e_eff)
        % Quasi-particles in XXZ-model are fermions
        F = -log(1+exp(-e_eff));
    end
    
    
    function h_i = getOneParticleEV(obj, t, x, rapid, i)
        switch i
        case 0 % eigenvalue of S_z operator
            h_i = repmat(obj.type_grid, obj.N, 1);
        case 1 
            h_i = obj.getBareMomentum(t, x, rapid,  obj.type_grid);
        case 2 % eigenvalue of H operator
            h_i = obj.getBareEnergyNoField(t, x, rapid,  obj.type_grid);
        otherwise 
            error(['Eigenvalue ' num2str(i) ' not implmented!'])
        end
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
        % Quasi-particles in XXZ-model are fermions
       theta = 1./(1 + exp(e_eff)); 
    end
    
    
    function [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid, type)        
        % Calculates velocities, outputs (1xNxM)
        de_dr   = obj.applyDressing(obj.calcEnergyRapidDeriv(t, x, rapid, type), theta, t);
        dp_dr   = obj.applyDressing(obj.calcMomentumRapidDeriv(t, x, rapid, type), theta, t);
        
        % Calculate acceleration from inhomogenous B-field. Note dmBdt
        % does not contribute as f = 0;
        a_eff_B = 0;
        if ~isempty(obj.couplings.dBdx)
            L_B     = repmat(type, obj.N, 1);
            a_eff_B = obj.couplings.dBdx(t,x).*obj.applyDressing(L_B, theta, t);
        end
        
        % Calculate acceleration from inhomogenous interaction
        a_eff_Delta = 0;
        if ~isempty(obj.couplings.dDeltadt) || ~isempty(obj.couplings.dDeltadx)
            % Calculate derivative of scattering phase with respect to
            % interaction c
            delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);            
            dT      = obj.calcScatteringDeltaDeriv(t, x, rapid, obj.rapid_grid, type, obj.type_grid);
            accKern = delta_k/(2*pi) * dT.*transpose(theta);
        end
        
        if ~isempty(obj.couplings.dDeltadt) % calc time deriv contribution
            f       = -obj.calcMomentumDeltaDeriv(t, x, rapid, type) + accKern*dp_dr;
            f_dr    = obj.applyDressing(f, theta, t);
            a_eff_Delta = a_eff_Delta + obj.couplings.dDeltadt(t,x).*f_dr;
        end

        if ~isempty(obj.couplings.dDeltadx) % calc space deriv contribution
            L       = -obj.calcEnergyDeltaDeriv(t, x, rapid, type) + accKern*de_dr;
            L_dr    = obj.applyDressing(L, theta, t);
            a_eff_Delta = a_eff_Delta + obj.couplings.dDeltadx(t,x).*L_dr;
        end
        
        v_eff   = de_dr./dp_dr;
        a_eff   = (a_eff_Delta + a_eff_B)./dp_dr;
    end
    
    
end % end private methods
    
end % end classdef