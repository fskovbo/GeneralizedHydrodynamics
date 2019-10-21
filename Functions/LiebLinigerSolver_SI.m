classdef LiebLinigerSolver_SI
    % Wraps around LLS such that one can pass and receive quantitites in SI
    % units.

    
properties (Access = public)
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
    
    LLS         = []; % LiebLinigerSolver

end % end private properties


methods (Access = public)
    
    %% Constructor
    function obj = LiebLinigerSolver_SI(omega_perp_si, x_grid, rapid_grid, couplings, Options)
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
        couplings   = obj.convert2TBA(couplings, 'couplings');
        
        if nargin < 5
            disp('Warning! Auto-derivative of couplings not supported for SI units! Disabling auto-deriv ...')
            Options.autoDerivCoup = false;
        elseif Options.autoDerivCoup
            disp('Warning! Auto-derivative of couplings not supported for SI units! Disabling auto-deriv ...')
            Options.autoDerivCoup = false;
        end
        
        obj.LLS     = LiebLinigerSolver(x_grid, rapid_grid, couplings, Options);
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
                quantity_tba{1,1} = @(t, x) quantity_si{1,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si; % mu is in units of energy
                quantity_tba{1,2} = @(t, x) quantity_si{1,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si; % c is in units of inverse length
                
                % time derivatives
                quantity_tba{2,1} = @(t, x) quantity_si{2,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.t_si; 
                quantity_tba{2,2} = @(t, x) quantity_si{2,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.t_si; 
                
                % space derivatives
                quantity_tba{3,1} = @(t, x) quantity_si{3,1}( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si*obj.Lg_si; 
                quantity_tba{3,2} = @(t, x) quantity_si{3,2}( t*obj.t_si, x*obj.Lg_si )*obj.Lg_si*obj.Lg_si; 
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
                % Couplings
                quantity_si{1,1} = @(t, x) quantity_tba{1,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si; % mu is in units of energy
                quantity_si{1,2} = @(t, x) quantity_tba{2,1}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si; % c is in units of inverse length
                
                % time derivatives
                quantity_si{2,1} = @(t, x) quantity_tba{2,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si/obj.t_si; 
                quantity_si{2,2} = @(t, x) quantity_tba{2,2}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si/obj.t_si; 
                
                % space derivatives
                quantity_si{3,1} = @(t, x) quantity_tba{3,1}( t/obj.t_si, x/obj.Lg_si )*obj.Eg_si/obj.Lg_si; 
                quantity_si{3,2} = @(t, x) quantity_tba{3,2}( t/obj.t_si, x/obj.Lg_si )/obj.Lg_si/obj.Lg_si; 
            otherwise
                error('Conversion of specified unit is not implemented!')
        end 
    end
    
    %% Wrapper functions
    function mu0_fit = fitAtomnumber(obj, T, V_ext, Natoms, setCouplingFlag)
        % Convert SI --> TBA
        T       = obj.convert2TBA(T, 'temperature');
        V_ext   = obj.convertFunction( V_ext, 1/obj.Eg_si, [ obj.t_si, obj.Lg_si ]);
%         V_ext   = @(t, x) V_ext( t*obj.t_si, x*obj.Lg_si )/obj.Eg_si;
        
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
%         u_t = obj.convert2SI(u_t, 'length');
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


methods (Access = private)
    
    function f_out = convertFunction(obj, f_in, f_scale, arg_scale)
        func_str = func2str(f_in); % convert function to string
        args_str = extractBetween(func_str, '@(', ')'); % get arguments 
        body_str = extractAfter( func_str, ['@(' args_str{1}, ')'] ); % get body of function
        args_list= split( args_str{1}, ',' ); % separate each argument in individual cells

        % Replace all arguments with their scaled version
        assert( length(args_list) == length(arg_scale) ) % #arguments equal to #argument scalings

        op = '[. * / ^ ( ) + -]'; % String with mathematical operations
        for i = 1:length(args_list)
            % Argument to be replaced
            expr = args_list{i};

            % Replacement is scaled argument
            repl = ['(' num2str(arg_scale(i)) '*'  expr ')'];

            % Replace all occurences of expr, where it is surrounded by
            % mathematical operators.
            % Ex: will replace the 'x' in '4*x', but not in 'exp'
            idx = regexp(body_str, [op expr op ]);
            for j = fliplr(1:length(idx))
                body_str = replaceBetween( body_str, idx(j)+1, idx(j)+length(expr), repl );
            end
        end

        % Add arguments at beginning and function scaling in the end
        func_str = ['@(' args_str{1}, ')' body_str '*' num2str(f_scale)];

        % Convert back to anonymous function
        f_out = str2func(func_str);
    end
    
end % end private methods
    
end % end classdef