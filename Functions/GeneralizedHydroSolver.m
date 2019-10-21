classdef GeneralizedHydroSolver < handle
    % Superclass for solving GHD dynamical equations.
    % Extends handle class --> properties can be changed with setter
    % methods while retaining GHD object.
    % This class contains abstract methods, which are model specific,
    % whereby they must be implemented in the corresponding subclass.
    % 
    % This class and any subclass of it must follow the following
    % convention for indeces:
    %   Index 1: Rapid index used only for convolutions.
    %   Index 2: Standard rapid index, dimension N.
    %   Index 3: Standard space index, dimension M.
    %   Index 4: Standard time index.
    %
    % Thus, quantites have dimensions (1xNxM), while kernels have (NxNxM)
    
properties (Access = public)
    % Grid lengths
    M               = []; % number of spatial grid-points
    N               = []; % number of rapidity grid-points
    Ntypes          = []; % number of quasi-particle types
    
    % Grids (all vectors)
    rapid_grid      = [];
    x_grid          = [];
    type_grid       = [];
    
    % Cell array with anonymous functions
    model           = []; % @(t,x,rapid,type) energy(first), momentum (second), and scattering (third)
    couplings       = []; % @(t,x) coupling #i
    couplingDerivs  = []; % @(t,x) coupling #i derivatives w.r.t. t (first row) and x (second row)
    modelCoupDerivs = []; % @(t,x,rapid,type) derivatives of energy(first row), momentum (second row), and scattering (third row) w.r.t. couplings #i'th column
    modelRapidDerivs= []; % @(t,x,rapid,type) derivatives of energy(first), momentum (second), and scattering (third ) w.r.t. rapidity
    
    % Optional parameters (default values specified here)
    homoEvol        = false;    % Flag for homogeneous evolution (no acceleration)
    tolerance       = 1e-6;     % Tolerance for TBA solution
    maxcount        = 100;      % Max interations for TBA solution
    stepOrder       = 2;        % Order of time step
    extrapFlag      = false;    % Flag for using extrapolation during propagation (useful for open systems)
    periodRapid     = false;    % Flag for periodic boundary conditions on rapidity
    simplifySteps   = 100;      % Number of simplification steps used for auto-deriv
    autoDerivModel  = true;     % Flag for using auto-deriv of model functions. If false, all corresponding methods must be overloaded!
    autoDerivCoup   = true; 	% Flag for using auto-deriv of couplings. If false, all couplings must be passed!

end % end protected properties


properties (Abstract, Access = protected)
    % All these properties must be specified in the specific implementation
    
    % Formulas of model formated as strings
    energy
    momentum
    scattering
    
    % Species of quasiparticle
    quasiSpecies
    
    % Names of couplings in cell (must match model formulas)
    couplingNames
    
end % end abstract properties


methods (Abstract, Access = protected)
    
    h_i     = getOneParticleEV(obj, t, x, rapid, i)

end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = GeneralizedHydroSolver(x_grid, rapid_grid, couplings, Ntypes, Options)        
        % Check format of input
        if ~isvector(x_grid) || ~isvector(rapid_grid)
            error('Input has wrong format!')
        end
        if ~iscell( couplings )
            error('couplings must be cell array of anonymous functions!')
        end
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.Ntypes      = Ntypes;
        
        % Reshape grids to right format
        obj.x_grid      = reshape(x_grid, 1, 1, 1, 1, obj.M); % 5th index is space
        obj.rapid_grid  = reshape(rapid_grid, obj.N, 1); % 1st index is rapidity
        obj.type_grid   = reshape( 1:Ntypes, 1, 1, Ntypes ); % Types are 3rd index
        
        % Copy fields of Options struct into class
        if nargin > 4
            fn = fieldnames(Options);
            for i = 1:length(fn)                      
                if isprop(obj, fn{i}) % only copy field if defined among properties
                    eval(['obj.',fn{i},' = Options.',fn{i},';']);
                end
            end
        end
        
        % Fill out cell arrays with anonymous functions
        obj.setCouplings(couplings);
        
    end
    
    
    function setCouplings(obj, couplings)
        % Sets couplings and computes all cell arrays with anonymous
        % founctions.
        if ~iscell( couplings )
            error('couplings must be cell array of anonymous functions!')
        end
        assert(size(couplings, 2) == length(obj.couplingNames))
        
        obj.couplings = couplings(1,:);
        
        
        if ~obj.autoDerivCoup
            % Derivatives of coupling are passed along with couplings
            obj.couplingDerivs = couplings([2,3],:);
        else
            % Calculate derivatives of couplings w.r.t. x and t
            obj.couplingDerivs = cell(2, length(couplings));

            for i = 1:length(couplings)           
                coup_str    = func2str( couplings{i} );
                coup_str    = erase(coup_str, ' '); % erase spaces
                coup_str    = erase(coup_str, '@(t,x)'); % erase function prefix
                [dcdt, z1]  = obj.takeDerivStr( coup_str, 't' );
                [dcdx, z2]  = obj.takeDerivStr( coup_str, 'x' );

                % If non-zero set derivs, otherwise leave empty
                if ~z1; obj.couplingDerivs{1,i} = obj.formatFunction(dcdt, false); end
                if ~z2; obj.couplingDerivs{2,i} = obj.formatFunction(dcdx, false); end
            end
        end
        
        
        % Format model functions to anonymous functions
        obj.model{1}    = obj.formatFunction(obj.energy, false);
        obj.model{2}    = obj.formatFunction(obj.momentum, false);
        obj.model{3}    = obj.formatFunction(obj.scattering, true);
        
        % Take derivatives of model w.r.t. rapidity and format to anonymous functions
        if obj.autoDerivModel
            dedr = obj.takeDerivStr( obj.energy, 'rapid' );
            dpdr = obj.takeDerivStr( obj.momentum, 'rapid' );
            dTdr = obj.takeDerivStr( obj.scattering, 'rapid' );

            obj.modelRapidDerivs{1} = obj.formatFunction(dedr, false);
            obj.modelRapidDerivs{2} = obj.formatFunction(dpdr, false);
            obj.modelRapidDerivs{3} = obj.formatFunction(dTdr, true);

            % Take derrivatives of model w.r.t. couplings and format to anonymous functions
            for i = 1:length(obj.couplingNames)
                dedc = obj.takeDerivStr( obj.energy, obj.couplingNames{i} );
                dpdc = obj.takeDerivStr( obj.momentum, obj.couplingNames{i} );
                dTdc = obj.takeDerivStr( obj.scattering, obj.couplingNames{i} );

                obj.modelCoupDerivs{1,i} = obj.formatFunction(dedc, false);
                obj.modelCoupDerivs{2,i} = obj.formatFunction(dpdc, false);
                obj.modelCoupDerivs{3,i} = obj.formatFunction(dTdc, true);
            end
        end
        
        
        % If all coupling derivatives are empty (equal to 0), set flag for
        % homogeneous evolution
        if all(all( cellfun(@isempty, obj.couplingDerivs) ))
            obj.homoEvol = true;
        else
            obj.homoEvol = false;
        end
    end
    
    
    function couplings = getCouplings(obj)
        if obj.autoDerivCoup
            couplings = obj.couplings;
        else
            % Include derivatives
            couplings = [obj.couplings ; obj.couplingDerivs];
        end
    end
    
    
    function [theta_t, u_t] = propagateTheta(obj, theta_init, t_array)
        % Propagates filling fraction according to GHD equation.
        % Simultaneously calculates characteristic function, u, used for
        % correlation functions.
        
        Nsteps          = length(t_array) - 1;
        
        theta_t         = cell(1, Nsteps+1);
        theta_t{1}      = theta_init;
        theta           = theta_init;
        
        u_t             = cell(1, Nsteps+1);
        u_init          = GHDtensor( repmat( obj.x_grid, obj.N, 1, obj.Ntypes, 1) );
        u_t{1}          = u_init;

        
        % Declare step as anonymous function according to stepOrder setting
        switch obj.stepOrder
            case 1
                % Set stepfunction
                stepFunc    = @(th_aux, th_prev, u_prev, t, dt) obj.performFirstOrderStep(th_aux, th_prev, u_prev, t, dt);                              
                
                % Calculate initial theta_aux (here theta_guess)
                theta_aux   = GHDtensor(obj.N, 1, obj.Ntypes, 1, obj.M);
                u           = u_init;
                
            case 2
                % Set stepfunction
                stepFunc    = @(th_aux, th_prev, u_prev, t, dt) obj.performSecondOrderStep(th_aux, th_prev, u_prev, t, dt);
                
                % Calculate initial th_aux (here theta_mid)
                dt          = t_array(2) - t_array(1);
                ddt         = dt/2/10;
                theta_aux   = theta_init;
                u           = u_init;
                theta_guess = GHDtensor(obj.N, 1, obj.Ntypes, 1, obj.M);
                
                % Calculate first theta_mid at dt/2 using first order step
                for i = 1:10
                    t         = (i-1)*ddt;
                    theta_aux = obj.performFirstOrderStep(theta_guess, theta_aux, u_init, t, ddt);
                end
                
            otherwise
                error(['Step order ' num2str(obj.stepOrder) ' not supported!'])
        end

        % initialize progress bar
        cpb = ConsoleProgressBar();                 % create instance
        initialize_progbar;                         % initialize progress bar with standard parameters
        fprintf('Time evolution progress:');
        cpb.start();   

        % Propagate theta using stepfunction
        for n = 1:Nsteps
            dt                  = t_array(n+1) - t_array(n);
            [theta, theta_aux, u] = stepFunc(theta_aux, theta, u, t_array(n), dt);
            theta_t{n+1}  = theta;
            u_t{n+1}      = u; 
            
            % show progress
            cpb_text = sprintf('%d/%d steps evolved', n, Nsteps);
            cpb.setValue(n/Nsteps);
            cpb.setText(cpb_text);
        end
        fprintf('\n')
    end
    
       
    function [rho, rhoS] = transform2rho(obj, theta, t_array)
        % Transforms from occupation number basis (theta) to
        % occupied density of states basis (rho).
        % NOTE: if couplings are time-dependent, t = 0 should be used for
        % units to stay fixed.
        
        Nsteps  = length(theta); % number of time steps
        rho     = cell(1,Nsteps);
        rhoS    = cell(1,Nsteps);
        for n = 1:Nsteps % Transform for each time step
            
            if Nsteps == 1
                theta_n = theta;
            else
                theta_n = theta{n};
            end
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            dp_dr   = obj.applyDressing(dp, theta_n, t);
            
            rhoS{n} = 1/(2*pi) * dp_dr;
            rho{n}  = theta_n.*rhoS{n}; 
        end
        
        if Nsteps == 1
            rho     = rho{1};
            rhoS    = rhoS{1};
        end
    end
    
    
    function [theta, rhoS] = transform2theta(obj, rho, t_array)
        % Transforms from occupied density of states basis (rho) to
        % occupation number basis (theta).
        % NOTE: if couplings are time-dependent, t = 0 should be used for
        % units to stay fixed.
        
        Nsteps  = length(rho); % number of time steps
        theta   = cell(1,Nsteps);
        rhoS    = cell(1,Nsteps);
        for n = 1:Nsteps % Transform for each time step
            if Nsteps == 1
                rho_n = rho;
            else
                rho_n = rho{n};
            end
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);
            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            kernel  = delta_k/(2*pi) * obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid); 
            
            rhoS{n} = dp/(2*pi) - kernel*rho_n;
            theta{n}= rho_n./rhoS{n};
        end
        
        if Nsteps == 1
            theta   = theta{1};
            rhoS    = rhoS{1};
        end
    end
    

    function [q, j] = calcCharges(obj, theta, c_idx, t_array)
        % Calculate charges, q_i, and associated currents, j_i, where i is
        % entry in vector c_idx.
        
        Nsteps  = length(theta); % number of time steps
        Ncharg  = length(c_idx); 
        q       = zeros(obj.M, Nsteps, Ncharg);
        j       = zeros(obj.M, Nsteps, Ncharg);
        delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);
    
        for i = 1:Nsteps
            if Nsteps == 1
                theta_i = theta;
            else
                theta_i = theta{i};
            end
            
            if nargin < 4
                t = 0;
            else
                t = t_array(i);
            end
            
            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);
            dE      = obj.calcEnergyRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.type_grid);

            for n = 1:Ncharg
                hn          = obj.getOneParticleEV( t, obj.x_grid, obj.rapid_grid, c_idx(n));               
                hn_dr       = obj.applyDressing(hn, theta_i, t);
                
                q(:,i,n)    = delta_k/(2*pi) * squeeze(sum( sum( double(dp.*theta_i.*hn_dr) , 3) , 1)); % (M,1)
                j(:,i,n)    = delta_k/(2*pi) * squeeze(sum( sum( double(dE.*theta_i.*hn_dr) , 3) , 1)); % (M,1)
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % Calculate initial filling function, given by the equilibrium TBA
        % for some temperature, T, and couplings.
        % If no couplings are passed to method, use already set couplings.
        
        if nargin == 3 % use input TBA_couplings
            couplings_old = obj.getCouplings(); % save old couplings
            obj.couplings = TBA_couplings(1,:);
%             obj.setCouplings(TBA_couplings);
        end
            
        e_eff = obj.calcEffectiveEnergy(T, 0, obj.x_grid, obj.rapid_grid);
        theta = obj.calcFillingFraction(e_eff);
        
        if nargin == 3
            obj.couplings = couplings_old(1,:);
%             obj.setCouplings(couplings_old); % return to old couplings
        end
    end
    

    function [CM, theta, C1P] = calcCorrelationMatrix(obj, T, init_couplings, tcorr_array, dt, c_idx, areCurrents)
        % Calculate two-point correlation matrix of two observables (O1,O2)
        %   CM(x,y) = <O1(x,t) O2(y,0)>
        % as well as one-point correlations, C1P, (i.e. averages) of said
        % observables. Here, times t are entries of tcorr_array
        %
        % Currently only charge-densities and their associated currents are
        % supported as observables. Thus,
        %   O1 = q_i (for areCurrents(1) = 0 and i = c_idx(1))
        %   O1 = j_i (for areCurrents(1) = 1 and i = c_idx(1))
        % and similarly for O2.
        %
        % NOTE: Initial state must be a TBA stationary state specified by
        % the parameters (T, init_couplings). Furthermore, any subsequent
        % evolution of the state must be homogeneous!        
        
        
        % Calculate initial state
        [theta_0, e_eff]    = obj.calcThermalState(T, init_couplings);
        [rho_0, rhoS_0]     = obj.transform2rho(theta_0, 0);
        
        % Calculate "acceleration" i.e. degree of inhomogeniety of state
        [~,~,~,~,dedx]  = gradient( double(e_eff), 0, 0, 0, 0 , permute(obj.x_grid, [5 1 2 3 4]) );
        a_eff0      = -dedx./(2*pi*rhoS_0);
        
        % Define all necessary functions for GHDcorrelations class
        GHDeqs.applyDressing            = @(Q, theta, t) obj.applyDressing(Q, theta, t);
        GHDeqs.calcScatteringRapidDeriv = @(rapid1, rapid2) obj.calcScatteringRapidDeriv(1, obj.x_grid, rapid1, rapid2, obj.type_grid, obj.type_grid);
        GHDeqs.getStatFactor            = @(theta) obj.getStatFactor(theta);
        
        % Initialize object for calculating the correlations
        GHDcorr = GHDcorrelations(obj.x_grid, obj.rapid_grid, theta_0, rho_0, rhoS_0, a_eff0, GHDeqs, obj.periodRapid);
        
        
        % Calculate state variables at correlation times.
        % NOTE: propagation MUST be homogeneous!
         
        % Create time grid for evolution
        t_array         = 0:dt:tcorr_array(end);
        t_array         = sort( [t_array, tcorr_array] ); % make sure tcorr_array is in t_array 
        t_array         = unique(t_array); % removes dublicates if tcorr_array was already in t_array
        [~, t_idx]      = ismember(tcorr_array, t_array); % get indices of tcorr_array entries
        
        % Time evolve theta_0
        [theta_t, u_t]  = obj.propagateTheta(theta_0, t_array);
        theta_t         = theta_t(t_idx); % get theta at t = 0 and t in tcorr_array
        u_t             = u_t(t_idx); 
        [rho_t, rhoS_t] = obj.transform2rho(theta_t, tcorr_array);
        
        
        % Calculate/prepare outputs
        [q0, j0]        = obj.calcCharges(theta_0, c_idx(2), 0);
        CM              = zeros(obj.M, obj.M, 2, length(tcorr_array));
        theta           = [ {theta_0}, theta_t];
        C1P             = zeros(obj.M, length(tcorr_array) + 1);
         
        for k = 1:length(tcorr_array)
            tic
            fprintf('Calculating correlation matrix (%d/%d) ... \n', k , length(tcorr_array));
            
            % Calculate g-functions, which will be acted on by the
            % correlation propagators.
            hi = obj.getOneParticleEV(tcorr_array(k), obj.x_grid, obj.rapid_grid, c_idx(1)); % should be (1xN)
            hj = obj.getOneParticleEV(0, obj.x_grid, obj.rapid_grid, c_idx(2));

            g_t = obj.applyDressing( hi, theta_t{k}, tcorr_array(k) ); % should be (1xNxM)
            g_0 = obj.applyDressing( hj, theta_0, 0 );
            
            [qt, jt] = obj.calcCharges(theta_t{k}, c_idx(1), tcorr_array(k));
            
            C1P(:,1)    = q0;
            C1P(:,1+k)  = qt;
            
            if areCurrents(1)
                g_t = g_t .* obj.calcEffectiveVelocities(theta_t{k}, tcorr_array(k), obj.x_grid, obj.rapid_grid, obj.type_grid);
                C1P(:,1+k)  = jt;
            end
            if areCurrents(2)
                g_0 = g_0 .* obj.calcEffectiveVelocities(theta_0, 0, obj.x_grid, obj.rapid_grid, obj.type_grid);
                C1P(:,1)    = j0;
            end
            
            
            % Calculate correlation matrix
            CM(:,:,:,k) = GHDcorr.calcCorrelationMatrix(g_0, g_t, theta_t{k}, rho_t{k}, rhoS_t{k}, u_t{k});
            
            toc
            fprintf('\n')
        end
    end
    
    
end % end public methods


methods (Access = protected)
    
    %% Functions for model functions and their derivatives
    % This section containsFunctions for calculating model parameters
    % (energy, momentum, scattering) along with their derivatives.
    % Note, these functions all call anonymous functions stored in cell
    % arrays, which are generated using the auto-differentiation feature.
    % If this feature does not produce good functions, one can manually
    % overload these functions in the specific implementation.
    function ebare = getBareEnergy(obj, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        ebare = obj.model{1}(t, x, rapid, type);
        
        if size(ebare, 1) == 1
            ebare = repmat(ebare, obj.N, 1);
        end
    end
    
    
    function pbare = getBareMomentum(obj, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        pbare = obj.model{2}(t, x, rapid, type);
        
        if size(pbare, 1) == 1
            pbare = repmat(pbare, obj.N, 1);
        end
    end
    
    
    function de = calcEnergyRapidDeriv(obj, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        de = obj.modelRapidDerivs{1}(t, x, rapid, type);
        
        if size(de, 1) == 1
            de = repmat(de, obj.N, 1);
        end
    end

    
    function dp = calcMomentumRapidDeriv(obj, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        dp = obj.modelRapidDerivs{2}(t, x, rapid, type);
        
        if size(dp, 1) == 1
            dp = repmat(dp, obj.N, 1);
        end
    end
    
    
    function dT = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2, type1, type2)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        dT = obj.modelRapidDerivs{3}(t, x, rapid1, rapid2, type1, type2);
        
        dT(isnan(dT)) = 0; % removes any NaN
        
        dT = GHDtensor(dT); % Converts to GHDtensor
    end
    
    
    function de = calcEnergyCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        de = obj.modelCoupDerivs{1,coupIdx}(t, x, rapid, type);
        
        if size(de, 1) == 1
            de = repmat(de, obj.N, 1);
        end
    end

    
    function dp = calcMomentumCouplingDeriv(obj, coupIdx, t, x, rapid, type)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        dp = obj.modelCoupDerivs{2,coupIdx}(t, x, rapid, type);
        
        if size(dp, 1) == 1
            dp = repmat(dp, obj.N, 1);
        end
    end
    
    
    function dT = calcScatteringCouplingDeriv(obj, coupIdx, t, x, rapid1, rapid2, type1, type2)
        % This function simply calls the appropriate anonnymous function
        % stored in one of the cell arrays, and makes sure that the format
        % adheres to convention.
        % Is coded as function for clear name and overload posibility.
        
        % Reshape input to right dimensions
        rapid1  = reshape(rapid1, length(rapid1), 1); % rapid1 is 1st index
        rapid2  = reshape(rapid2, 1, length(rapid2)); % rapid2 is 2nd index
        type1   = reshape(type1, 1, 1, length(type1)); % type1 is 3rd index
        type2   = reshape(type2, 1, 1, 1, length(type2)); % type2 is 4th index
        
        dT = obj.modelCoupDerivs{3,coupIdx}(t, x, rapid1, rapid2, type1, type2);
        
        dT(isnan(dT)) = 0; % removes any NaN
        
        dT = GHDtensor(dT); % Converts to GHDtensor
    end
    
    
    function f = getStatFactor(obj, theta)
        switch obj.quasiSpecies
        case 'fermion'
            f = 1 - theta;
        case 'boson'
            f = 1 + theta;
        case 'classical'
            f = 1;
        case 'radiative'
            f = theta;
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end
    end

    
    function F = getFreeEnergy(obj, e_eff)
        switch obj.quasiSpecies
        case 'fermion'
            F = -log( 1 + exp(-e_eff));
        case 'boson'
            F = log( 1 - exp(-e_eff));
        case 'classical'
            F = -exp(-e_eff);
        case 'radiative'
            F = log( e_eff );
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end 
    end
    
    
    %% Standard GHD equations used for the public functions
    function Q_dr = applyDressing(obj, Q, theta, t)
        % Applies dressing to quantity Q. Q must have first index of
        % dimension N!
        % NOTE: in correlation paper a distinction is made between vector
        % and scalar field. By these conventions T (diff scattering) is
        % symmetric, whereby scalar and vector fields are treated equally.
        
        
        if ~isa(Q, 'GHDtensor')
            Q = GHDtensor(Q);
        end
        
        % Calculate dressing operator
        delta_k = obj.rapid_grid(2) - obj.rapid_grid(1); 
        kernel  = delta_k/2/pi*obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid);
        
        I_rapid = eye(obj.N);
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        identity= I_rapid.*I_type;

        U       = identity + kernel.*transpose(theta); 
        
        % We now have the equation Q = U*Q_dr. Therefore we solve for Q_dr
        % using the '\' operation.
        Q_dr     = U\Q;
           
    end
    
    
    function theta = calcFillingFraction(obj, e_eff)
        switch obj.quasiSpecies
        case 'fermion'
            theta = 1./( exp(e_eff) + 1);
        case 'boson'
            theta = 1./( exp(e_eff) - 1);
        case 'classical'
            theta = exp(-e_eff);
        case 'radiative'
            theta = 1./exp(e_eff);
        otherwise
            error(['Quasi-particle species ' obj.quasiSpecies ' is not implemented! Check spelling and cases!'])
        end
    end
    
    
    function [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid, type)        
        % Calculates velocities
        de_dr   = obj.applyDressing(obj.calcEnergyRapidDeriv(t, x, rapid, type), theta, t);
        dp_dr   = obj.applyDressing(obj.calcMomentumRapidDeriv(t, x, rapid, type), theta, t);
        
        v_eff   = de_dr./dp_dr;
         
        % Calculate acceleration from inhomogenous couplings
        a_eff   = 0;
        
        if obj.homoEvol % if homogeneous couplings, acceleration = 0
            a_eff = GHDtensor( zeros(size( v_eff )) );
            return
        end
        
        % Calculate contribution for each coupling
        for coupIdx = 1:length(obj.couplings)
            delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);            
            dT      = obj.calcScatteringCouplingDeriv(coupIdx, t, x, rapid, obj.rapid_grid, type, obj.type_grid);
            accKern = delta_k/(2*pi) * dT.*transpose(theta);
            
            % if time deriv of coupling exist compute f
            if ~isempty(obj.couplingDerivs{1,coupIdx}) 
                f       = -obj.calcMomentumCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*dp_dr;
                f_dr    = obj.applyDressing(f, theta, t);
                a_eff   = a_eff + obj.couplingDerivs{1,coupIdx}(t,x).*f_dr;
            end
            
            % if spacial deriv of coupling exist compute Lambda
            if ~isempty(obj.couplingDerivs{2,coupIdx}) 
                L       = -obj.calcEnergyCouplingDeriv(coupIdx, t, x, rapid, type) + accKern*de_dr;
                L_dr    = obj.applyDressing(L, theta, t);
                a_eff   = a_eff + obj.couplingDerivs{2,coupIdx}(t,x).*L_dr;
            end
        end
        
     
        a_eff   = a_eff./dp_dr;
    end
    
    
    function e_eff = calcEffectiveEnergy(obj, T, t, x, rapid)
        % Calculates the dressed energy per particle epsilon(k), from which the
        % preasure and later the filling factor theta can be derived.
        % This is achieved by iteratively solving the TBA equation.
        dk          = obj.rapid_grid(2) - obj.rapid_grid(1);
        ebare       = obj.getBareEnergy(t, x, obj.rapid_grid, obj.type_grid); 
        kernel      = dk/(2*pi)*obj.calcScatteringRapidDeriv(t, x, obj.rapid_grid, obj.rapid_grid, obj.type_grid, obj.type_grid );
        
        e_eff       = GHDtensor(obj.N, 1, obj.Ntypes, 1, obj.M);
        e_eff_old   = GHDtensor(obj.N, 1, obj.Ntypes, 1, obj.M);
        error_rel   = 1;
        count       = 0;
        
        % Solve TBA eq. for epsilon(k) by iteration:
        % Using epsilonk_old, update epsilon_k until convergence is reached 
        % (difference between epsilonk_old and epsilon_k becomes less than tol)
        while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
            
            % calculate epsilon(k) from integral equation using epsilonk_old
            % i.e. update epsilon^[n] via epsilon^[n-1]            
            e_eff       = ebare./T - kernel*obj.getFreeEnergy(e_eff_old);
            
            % calculate error
            v1          = flatten(e_eff);
            v2          = flatten(e_eff_old);
            
            sumeff      = sum( v1.^2 ,1);            
            error_rel   = squeeze(sum( (v1 - v2).^2, 1)./sumeff);
            e_eff_old   = e_eff;

            count       = count+1;
        end
    end
    
    
    %% Time steppers and associated supporting functions
    function [theta_next, theta_tmp, u_next] = performFirstOrderStep(obj, theta, theta_prev, u_prev, t, dt)            
        % Use single Euler step
        [v_eff, a_eff]  = obj.calcEffectiveVelocities(theta_prev, t, obj.x_grid, obj.rapid_grid, obj.type_grid); % should be (1xNxM)
            
        x_back          = obj.x_grid - dt*v_eff;
        r_back          = obj.rapid_grid - dt*a_eff;

        % Use interpolation to find theta_prev at x_back, r_back and
        % assign values to theta.
        theta_next      = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
        u_next          = obj.interpPhaseSpace(u_prev, r_back, x_back, true ); % always extrapolate u
        theta_tmp       = zeros(obj.N,1,obj.Ntypes,1,obj.M); % Unused outout for similarity to higher order steps       
    end
    
    
    function [theta_next, theta_mid, u_next] = performSecondOrderStep(obj, theta_mid, theta_prev, u_prev, t, dt)
        % Calculate theta^[n+1] using the filling (theta) at the midpoint
        % between theta^[n] and theta^[n+1].
        [theta_next, u_next] = step2(obj, theta_mid, theta_prev, u_prev, t, dt);
        
        % Perform another step with newly calculated theta as midpoint, in
        % order to calculate the midpoint filling for the next step.
        % Doesnt output u_mid, as it is not needed.
        theta_mid  = step2(obj, theta_next, theta_mid, u_prev, t+dt/2, dt);
        
        
        function [theta_next, u_next] = step2(obj, theta_mid, theta_prev, u_prev, t, dt)
            % Estimate x' and rapid' using midpoint filling
            [v_eff, a_eff] = obj.calcEffectiveVelocities(theta_mid, t+dt/2, obj.x_grid, obj.rapid_grid, obj.type_grid);

            x_mid   = obj.x_grid - 0.5*dt*v_eff; % (1xNxM)
            r_mid   = obj.rapid_grid - 0.5*dt*a_eff; % (1xNxM)

            % Interpolate v_eff and a_eff to midpoint coordinates x' and rapid' 
            v_mid   = obj.interpPhaseSpace( v_eff, r_mid, x_mid, true ); % always extrapolate v_eff
            a_mid   = obj.interpPhaseSpace( a_eff, r_mid, x_mid, true );
            
            % From midpoint velocity calculate fullstep x and rapid translation
            x_back  = obj.x_grid - dt*v_mid;
            r_back  = obj.rapid_grid - dt*a_mid;

            % Use interpolation to find theta_prev at x_back, r_back and
            % assign values to theta_next.
            theta_next = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
            u_next     = obj.interpPhaseSpace(u_prev, r_back, x_back, true); % always extrapolate u
        end % end nested function
    end
    
    
    function tensor_int = interpPhaseSpace(obj, tensor_grid, rapid_int, x_int, extrapFlag)
        % This function exists because MATLAB has different syntax between
        % interp1 and interp2, and I want to choose whether i extrapolate
        % or not with a simple TRUE/FAlSE argument.
        % ASSUME function_grid is on x_grid and rapid_grid.
        % Returns function_int with same dimensions as input function_grid
        
        % rapid_int and x_int should be (N,1,Nt,1,M)
        
        % Cast to matrix form
        x_int       = double(x_int);
        rapid_int   = double(rapid_int);
        mat_grid    = double(tensor_grid); % should be (N,1,Nt,1,M)
        
        % Need spacial dimension as first index in order to use (:) linearization
        x_int       = permute(x_int, [5 1 3 4 2]); % (M,N,Nt,1,1)
        rapid_int   = permute(rapid_int, [5 1 3 4 2]);
        
        x_g         = permute(obj.x_grid, [5 1 3 4 2]); % (M,N,Nt,1,1)
        rapid_g     = permute(obj.rapid_grid, [5 1 3 4 2]);
        
        % Enforce periodic boundary conditions
        if obj.periodRapid 
            rapid_int = mod(rapid_int + obj.rapid_grid(1), obj.rapid_grid(end)-obj.rapid_grid(1)) + obj.rapid_grid(1);
        end
        
        % Get matrix representation of GHDtensor and pemute spacial index
        % to first.
        mat_grid    = permute(mat_grid, [5 1 3 4 2]);
        mat_int     = zeros(obj.M, obj.N, obj.Ntypes);
        
        for i = 1:obj.Ntypes
            rapid_i = rapid_int(:,:,i);
            x_i     = x_int(:,:,i);
            mat_g   = mat_grid(:,:,i);   
            
            if extrapFlag
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), 'spline');
            else
                % Set all extrapolation values to zero!
                mat_tmp = interp2( rapid_g, x_g, mat_g, rapid_i(:), x_i(:), 'spline', 0);
            end
           
            mat_tmp = reshape(mat_tmp, obj.M, obj.N);
            mat_int(:,:,i) = mat_tmp;
        end
        
        % Add dummy indices and reshape back to original indices
        mat_int = permute(mat_int, [2 5 3 4 1] );
        
        tensor_int = GHDtensor(mat_int);
    end
  
    
    %% Functions for auto-derivatives
    function [deriv_str, isZero] = takeDerivStr(obj, func_str, var_str)
        % Takes a function formated as a string and output its derivative
        % (with respect to the string variable var_str) as a string

        % Convert to symbolics
        func_sym    = str2sym( func_str );
        var_sym     = str2sym( var_str );
        deriv_sym   = diff( func_sym, var_sym );

        % replace pi with placeholder symbolic, as simplify always treats pi as
        % a numeric value
        deriv_str   = char( deriv_sym );
        deriv_str   = regexprep(deriv_str, 'pi', 'placeholder');
        deriv_sym   = str2sym( deriv_str );

        deriv_sym   = simplify( deriv_sym , 'Steps', obj.simplifySteps);

        % If derivative is zero, return empty 
        if isequal(deriv_sym, sym(0))
            isZero = true;
        else
            isZero = false;
        end

        % Format back to string
        deriv_str   = char( deriv_sym );
        deriv_str   = regexprep(deriv_str, 'placeholder', 'pi');
    end
    
    
    function func = formatFunction(obj, func_str, isKernel)
        % Takes a function formated as a string and outputs an anonymous function 
          
        % Replace operations with elementwise ones
        func_str = regexprep(func_str, '*', '.*');
        func_str = regexprep(func_str, '/', './');
        func_str = regexprep(func_str, '\^', '\.^');
        
        % Replace all couplings with their expressions
        for i = 1:length(obj.couplingNames)
            if obj.autoDerivCoup
                % Replace couplings with their expressions
                coup_str = func2str( obj.couplings{i} );
                coup_str = erase(coup_str, ' '); % erase spaces
                coup_str = erase(coup_str, '@(t,x)'); % erase function prefix
                coup_str = [ '(' coup_str ')' ]; % brace expression
                func_str = regexprep(func_str, obj.couplingNames{i}, coup_str);
            else
                % Replace couplings with link to their anonymous function
                op = '[* / ^ ( ) + - \s]';
                repl = ['(obj.couplings{' num2str(i) '}(t,x))'];

                idx1 = regexp(func_str, [op obj.couplingNames{i} op ]);
                idx2 = regexp(func_str, [op obj.couplingNames{i} ]);
                idx3 = regexp(func_str, [ obj.couplingNames{i} op ]);
                idx  = unique( [ idx1 idx2 (idx3-1) ] );

                for j = fliplr(1:length(idx))
                    func_str = replaceBetween( func_str, idx(j)+1, idx(j)+length(obj.couplingNames{i}), repl );
                end
            end
        end
        
        if isKernel
            func_str = regexprep(func_str, 'rapid', '(rapid1-rapid2)');
            
            % Add parameter prefix to string 
            func_str = ['@(t, x, rapid1, rapid2, type1, type2) ' func_str];
        else
            % Add parameter prefix to string 
            func_str = ['@(t, x, rapid, type) ' func_str];
        end

        % Convert to anonymous function 
        if obj.autoDerivCoup
            func     = str2func(func_str);
        else
            % str2func doesnt work when string includes links to other
            % anonymous functions...
            func = eval(func_str);
        end
    end
    
end % end protected methods

end % end classdef