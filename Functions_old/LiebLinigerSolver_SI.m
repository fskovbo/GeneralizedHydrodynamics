classdef LiebLinigerSolver_SI
    % Wraps around LLS such that one can pass and receive quantitites in SI
    % units.

    
properties (Access = private)
    % Physical constants
    m_si        = 87*1.6605402e-27; % Rb-87 mass
    hbar_si     = 1.054571726e-34;
    kB_si       = 1.38065e-23;
    as_si       = 5.2e-9; % Rb-87 scattering length
    
    % TBA unit scales     
    Eg_si       = []; % energy scale
    Lg_si       = []; % length scale
    t_si        = []; % time scale
    T_si        = []; % temperature scale
    P_si        = []; % momentum scale
    
    LLcouplings = {'mu', 'dmudt', 'dmudx', 'c', 'dcdt', 'dcdx'};
    
    LLS         = []; % LiebLinigerSolver

end % end private properties


methods (Access = public)
    
    %% Constructor
    function obj = LiebLinigerSolver_SI(omega_perp_si, x_grid, rapid_grid, couplings, stepOrder, extrapFlag)
        % Calculate unit scales
        g1D_si      = 2*obj.hbar_si*omega_perp_si*obj.as_si;

        obj.Eg_si   = 0.5*obj.m_si*g1D_si^2 /obj.hbar_si^2; 
        obj.Lg_si   = obj.hbar_si^2 /(g1D_si*obj.m_si); 
        obj.t_si    = obj.hbar_si/obj.Eg_si; 
        obj.T_si    = obj.Eg_si/obj.kB_si;
        obj.P_si    = obj.hbar_si/obj.Lg_si;
        
        if isnan(omega_perp_si) % only for testing
            obj.Eg_si   = 1; 
            obj.Lg_si   = 1; 
            obj.t_si    = 1; 
            obj.T_si    = 1;
            obj.P_si    = 1;   
        end
        
        % Instantiate LiebLinigerSolver
        x_grid      = obj.convert2TBA(x_grid, 'length');
        rapid_grid  = obj.convert2TBA(rapid_grid, 'rapidity');
        couplings   = obj.convert2TBA(couplings, 'couplings'); % REWRITE THIS!
        
        if nargin < 5
            % Propagate theta using stepper of order X.
            stepOrder = 2;
        end
        if nargin < 6
            % Extrapolate theta when propagating (necessary for homogeneous states)
            extrapFlag = false;
        end
        
        obj.LLS     = LiebLinigerSolver(x_grid, rapid_grid, couplings, stepOrder, extrapFlag);
    end
    
    
    %% Unit-convertion functions
    function quantity_tba = convert2TBA(obj, quantity_si, unit)
        switch unit
            case 'energy'
                quantity_tba = quantity_si/obj.Eg_si;
            case 'rapidity'
                quantity_tba = quantity_si*obj.Lg_si;
            case 'momentum'
                quantity_tba = quantity_si/obj.P_si;
            case 'time'
                quantity_tba = quantity_si/obj.t_si;
            case 'length'
                quantity_tba = quantity_si/obj.Lg_si;
            case 'temperature'
                quantity_tba = quantity_si/obj.T_si;
            case 'couplings'
                % Anonymous functions are in SI units and take SI
                % arguments. Thus, convert arguments to SI and output to
                % TBA units.
                quantity_tba.mu     = @(t, x) quantity_si.mu( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si;
                quantity_tba.dmudx  = @(t, x) quantity_si.dmudx( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.Lg_si;
                quantity_tba.dmudt  = @(t, x) quantity_si.dmudt( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.t_si;
                quantity_tba.c      = @(t, x) quantity_si.c( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si; % c is in units of inverse length
                quantity_tba.dcdx   = @(t, x) quantity_si.dcdx( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.Lg_si;
                quantity_tba.dcdt   = @(t, x) quantity_si.dcdt( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.t_si;
                
                % Make sure empty couplings stay empty!  
                fn = fieldnames(quantity_si);
                for i = 1:length(fn)                      
                    if isempty( quantity_si.(fn{i}) )
                        quantity_tba.(fn{i}) = [];
                    end
                end
            otherwise
                error('Conversion of specified unit is not implemented!')
        end 
    end
    
    
    function quantity_si = convert2SI(obj, quantity_tba, unit)
        switch unit
            case 'energy'
                quantity_si = quantity_tba*obj.Eg_si;
            case 'momentum'
                quantity_si = quantity_tba*obj.P_si;
            case 'rapidity'
                quantity_si = quantity_tba/obj.Lg_si;
            case 'time'
                quantity_si = quantity_tba*obj.t_si;
            case 'length'
                quantity_si = quantity_tba*obj.Lg_si;
            case 'temperature'
                quantity_si = quantity_tba*obj.T_si;
            case 'couplings'
                % Anonymous functions are in SI units and take SI
                % arguments. Thus, convert arguments to SI and output to
                % TBA units.
                quantity_si.mu     = @(t, x) quantity_tba.mu( t*obj.t_si, x*obj.Lg_si )*obj.Eg_si;
                quantity_si.dmudx  = @(t, x) quantity_tba.dmudx( t*obj.t_si, x*obj.Lg_si )*obj.Eg_si/obj.Lg_si;
                quantity_si.dmudt  = @(t, x) quantity_tba.dmudt( t*obj.t_si, x*obj.Lg_si )*obj.Eg_si/obj.t_si;
                quantity_si.c      = @(t, x) quantity_tba.c( t*obj.t_si, x*obj.Lg_si )/obj.Lg_si; % c is in units of inverse length
                quantity_si.dcdx   = @(t, x) quantity_tba.dcdx( t*obj.t_si, x*obj.Lg_si )/obj.Lg_si/obj.Lg_si;
                quantity_si.dcdt   = @(t, x) quantity_tba.dcdt( t*obj.t_si, x*obj.Lg_si )/obj.Lg_si/obj.t_si;
                
                % Make sure empty couplings stay empty!  
                fn = fieldnames(quantity_si);
                for i = 1:length(fn)                      
                    if isempty( quantity_si.(fn{i}) )
                        quantity_tba.(fn{i}) = [];
                    end
                end
            otherwise
                error('Conversion of specified unit is not implemented!')
        end 
    end
    
    %% Wrapper functions
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, setCouplingFlag)
        % Convert SI --> TBA
        T       = obj.convert2TBA(T, 'temperature');
        V_ext   = @(t, x) V_ext( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si;
        
        % Run LLS function
        mu0_fit = obj.LLS.fitAtomnumber(T, V_ext, Natoms, setCouplingFlag);
        
        % Convert TBA --> SI
        mu0_fit = obj.convert2SI(mu0_fit, 'energy');
    end
    
    
    function setCouplings(obj, couplings)
        % Convert SI --> TBA
        couplings = obj.convert2TBA(couplings, 'couplings');
        
        obj.LLS.setCouplings(couplings);
    end
    
    
    function couplings = getCouplings(obj)
        couplings = obj.LLS.getCouplings();
        
        % Convert TBA --> SI
        couplings = obj.convert2SI(couplings, 'couplings');
    end
    
    
    function [theta_t, u_t] = propagateTheta(obj, theta_init, t_array)
        % Convert SI --> TBA
        t_array = obj.convert2TBA(t_array, 'time');
        
        % Run LLS function
        [theta_t, u_t] = obj.LLS.propagateTheta(theta_init, t_array);
        
        % Convert TBA --> SI
        u_t = obj.convert2SI(u_t, 'length');
    end
    
    
    function [rho, rhoS] = transform2rho(obj, theta, t_array)
        if nargin == 3       
            % Convert SI --> TBA
            t_array = obj.convert2TBA(t_array, 'time');

            % Run LLS function
            [rho, rhoS] = obj.LLS.transform2rho(theta, t_array);
        else
            % Run LLS function
            [rho, rhoS] = obj.LLS.transform2rho(theta);
        end
    end
    
    
    function [theta, rhoS] = transform2theta(obj, rho, t_array)
        % Convert SI --> TBA (for inverse quantities, use SI rather than TBA)
        rho     = obj.convert2SI( rho, 'length'); % convert 'per length' 
        rho     = obj.convert2SI( rho, 'rapidity'); % convert 'per rapidity' 
        
        if nargin == 3       
            % Convert SI --> TBA
            t_array = obj.convert2TBA(t_array, 'time');

            % Run LLS function
            [theta, rhoS] = obj.LLS.transform2theta(rho, t_array);
        else
            % Run LLS function
            [theta, rhoS] = obj.LLS.transform2theta(rho);
        end
    end
    
    
    function [q, j] = calcCharges(obj, theta, c_idx, t_array)
        % Convert SI --> TBA
        t_array = obj.convert2TBA(t_array, 'time');
        
        % Run LLS function
        [q, j] = obj.LLS.calcCharges(theta, c_idx, t_array);
        
        % Convert TBA --> SI
        % NOTE: Doesn't convert currents
        for i = length(c_idx)
            switch c_idx(i)
            case 0 % atomic density
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            case 1 % momentum density
                q(:,:,i) = obj.convert2SI(q(:,:,i), 'momentum'); 
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            case 2 %energy density
                q(:,:,i) = obj.convert2SI(q(:,:,i), 'energy'); 
                q(:,:,i) = obj.convert2TBA(q(:,:,i), 'length'); % convert 'per length'
                
            otherwise
                disp(['Warning: No known unit of charge nr. ' num2str(c_idx(i))])
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % Convert SI --> TBA
        T = obj.convert2TBA(T, 'temperature');
        
        if nargin == 3
            TBA_couplings   = obj.convert2TBA(TBA_couplings, 'couplings');
            [theta, e_eff]  = obj.LLS.calcThermalState(T, TBA_couplings);
        else
            [theta, e_eff]  = obj.LLS.calcThermalState(T);
        end
    end
    
    
end % end public methods
    
end % end classdef