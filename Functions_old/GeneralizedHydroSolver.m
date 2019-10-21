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
    
properties (Access = protected)
    M           = [];
    N           = [];
    rapid_grid  = [];
    x_grid      = [];
    
    couplings   = [];
    
    tolerance   = 1e-6;
    maxcount    = 100;
    stepOrder   = 0;
    extrapFlag  = false;

end % end protected properties


methods (Abstract, Access = protected)

    ebare   = getBareEnergy(obj, t, x, rapid);
    pbare   = getBareMomentum(obj, t, x, rapid);
    
    dT      = calcScatteringRapidDeriv(obj, t, x, rapid1, rapid2)
    de      = calcEnergyRapidDeriv(obj, t, x, rapid)
    dp      = calcMomentumRapidDeriv(obj, t, x, rapid)
    
    f       = getStatFactor(obj, theta)
    h_i     = getOneParticleEV(obj, i, rapid)
    theta   = calcFillingFraction(obj, e_eff)

    [v_eff, a_eff] = calcEffectiveVelocities(obj, theta, t, x, rapid)
  
end % end protected abstract methods


methods (Access = public)
    
    % Superclass constructor
    function obj = GeneralizedHydroSolver(x_grid, rapid_grid, couplings, stepOrder, extrapFlag)        
        % Check dimensionality of input
        if ~isvector(x_grid) || ~isvector(rapid_grid)
            error('Input has wrong format!')
        end
        
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.stepOrder   = stepOrder;
        obj.extrapFlag  = extrapFlag;
        
        % Reshape grids to right format
        obj.x_grid      = reshape(x_grid, 1, 1, obj.M); % 3rd index is space
        obj.rapid_grid  = reshape(rapid_grid, 1, obj.N); % 2nd index is rapidity
        
        obj.setCouplings(couplings);
    end
    
    
    function setCouplings(obj, couplings)
        obj.couplings = couplings; 
    end
    
    
    function couplings = getCouplings(obj)
        couplings = obj.couplings;
    end
    
    
    function [theta_t, u_t] = propagateTheta(obj, theta_init, t_array)
        % Propagates filling fraction according to GHD equation.
        % Simultaneously calculates characteristic function, u, used for
        % correlation functions.
        
        % check dimensionality
        if any(size(theta_init) ~= [1, obj.N, obj.M])
            error('theta_init does not have the right format!')
        end
        
        Nsteps          = length(t_array) - 1;
        
        theta_t         = zeros( 1, obj.N, obj.M, Nsteps+1 );
        theta_t(:,:,:,1)= theta_init;
        theta           = theta_init;
        
        u_t             = zeros( 1, obj.N, obj.M, Nsteps+1 );
        u_init          = repmat( obj.x_grid, 1, obj.N);
        u_t(:,:,:,1)    = u_init;
        
        % Declare step as anonymous function according to stepOrder setting
        switch obj.stepOrder
            case 1
                % Set stepfunction
                stepFunc    = @(th_aux, th_prev, u_prev, t, dt) obj.performFirstOrderStep(th_aux, th_prev, u_prev, t, dt);                              
                
                % Calculate initial theta_aux (here theta_guess)
                theta_aux   = zeros(1, obj.N, obj.M);
                u           = u_init;
                
            case 2
                % Set stepfunction
                stepFunc    = @(th_aux, th_prev, u_prev, t, dt) obj.performSecondOrderStep(th_aux, th_prev, u_prev, t, dt);
                
                % Calculate initial th_aux (here theta_mid)
                dt          = t_array(2) - t_array(1);
                ddt         = dt/2/10;
                theta_aux   = theta_init;
                u           = u_init;
                theta_guess = zeros(1, obj.N, obj.M);
                
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
            theta_t(:,:,:,n+1)  = theta;
            u_t(:,:,:,n+1)      = u; 
            
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
        
        [i1, i2, i3 , i4] = size(theta);        

        % check dimensionality
        if i2 ~= obj.N || i1 ~= 1 || i3 ~= obj.M
            error('theta does not have right format!')
        end

        rho     = zeros(size(theta));
        rhoS    = zeros(size(theta));
        for n = 1:i4 % Transform for each time step
            theta_n = theta(:,:,:,n);
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid);
            dp_dr   = obj.applyDressing(dp, theta_n, t);
            
            rhoS(:,:,:,n) = dp_dr/2/pi;
            rho(:,:,:,n) = theta_n.*rhoS(:,:,:,n); 
        end
    end
    
    
    function [theta, rhoS] = transform2theta(obj, rho, t_array)
        % Transforms from occupied density of states basis (rho) to
        % occupation number basis (theta).
        % NOTE: if couplings are time-dependent, t = 0 should be used for
        % units to stay fixed.
        
        [i1, i2, i3 , i4] = size(rho);        

        % check dimensionality
        if i2 ~= obj.N || i1 ~= 1 || i3 ~= obj.M
            error('theta does not have right format!')
        end

        theta   = zeros(size(rho));
        rhoS    = zeros(size(rho));
        for n = 1:i4 % Transform for each time step
            rho_n = rho(:,:,:,n);
            
            if nargin < 3
                t = 0;
            else
                t = t_array(n);
            end

            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid);
            dT      = obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid);
            delta_k = obj.rapid_grid(2) - obj.rapid_grid(1); 
            
            rhoS(:,:,:,n)  = (1/(2*pi)) * ( dp - delta_k*sum(dT.*permute(rho_n,[2 1 3]) , 1));
            theta(:,:,:,n) = rho_n./rhoS(:,:,:,n);
        end
    end
    
    
    function [q, j] = calcCharges(obj, theta, c_idx, t_array)
        % Calculate charges, q_i, and associated currents, j_i, where i is
        % entry in vector c_idx.
        
        Nsteps  = size(theta,4); % number of time steps
        Ncharg  = length(c_idx); 
        q       = zeros(obj.M, Nsteps, Ncharg);
        j       = zeros(obj.M, Nsteps, Ncharg);
        dr      = obj.rapid_grid(2) - obj.rapid_grid(1);
    
        for i = 1:Nsteps 
            theta_i = theta(:,:,:,i);
            
            if nargin < 4
                t = 0;
            else
                t = t_array(i);
            end
            
            dp      = obj.calcMomentumRapidDeriv(t, obj.x_grid, obj.rapid_grid);
            dE      = obj.calcEnergyRapidDeriv(t, obj.x_grid, obj.rapid_grid);

            for n = 1:Ncharg
                hn          = obj.getOneParticleEV( c_idx(n), obj.rapid_grid);               
                hn_dr       = obj.applyDressing(hn, theta_i, t);
                
                q(:,i,n)    = dr/2/pi * reshape(reshape(permute(theta_i.*hn_dr,[2,1,3]),obj.N,[]).'*dp',obj.M,[]); % (M,1)
                j(:,i,n)    = dr/2/pi * reshape(reshape(permute(theta_i.*hn_dr,[2,1,3]),obj.N,[]).'*dE',obj.M,[]); % (M,1)
            end
        end
    end
    
    
    function [theta, e_eff] = calcThermalState(obj, T, TBA_couplings)
        % Calculate initial filling function, given by the equilibrium TBA
        % for some temperature, T, and couplings.
        % If no couplings are passed to method, use already set couplings.
        
        if nargin == 3 % use input TBA_couplings
            couplings_old = obj.getCouplings(); % save old couplings
            obj.setCouplings(TBA_couplings);
        end
            
        e_eff = obj.calcEffectiveEnergy(T, 0, obj.x_grid, obj.rapid_grid);
        theta = obj.calcFillingFraction(e_eff);
        
        if nargin == 3
            obj.setCouplings(couplings_old); % return to old couplings
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
        [~,~,dedx]          = gradient( e_eff, 0, 0, permute(obj.x_grid, [3 2 1]) );
        a_eff0              = -dedx./(2*pi*rhoS_0);
        
        % Define all necessary functions for GHDcorrelations class
        GHDeqs.applyDressing            = @(Q, theta, t) obj.applyDressing(Q, theta, t);
        GHDeqs.calcScatteringRapidDeriv = @(rapid1, rapid2) obj.calcScatteringRapidDeriv(1, obj.x_grid, rapid1, rapid2);
        GHDeqs.getStatFactor            = @(theta) obj.getStatFactor(theta);
        
        % Initialize object for calculating the correlations
        GHDcorr = GHDcorrelations(obj.x_grid, obj.rapid_grid, theta_0, rho_0, rhoS_0, a_eff0, GHDeqs);
        
        
        % Calculate state variables at correlation times.
        % NOTE: propagation MUST be homogeneous!
         
        % Create time grid for evolution
        t_array         = 0:dt:tcorr_array(end);
        t_array         = sort( [t_array, tcorr_array] ); % make sure tcorr_array is in t_array 
        t_array         = unique(t_array); % removes dublicates if tcorr_array was already in t_array
        [~, t_idx]      = ismember(tcorr_array, t_array); % get indices of tcorr_array entries
        
        % Time evolve theta_0
        [theta_t, u_t]  = obj.propagateTheta(theta_0, t_array);
        theta_t         = theta_t(:,:,:,t_idx); % get theta at t = 0 and t in tcorr_array
        u_t             = u_t(:,:,:,t_idx); 
        [rho_t, rhoS_t] = obj.transform2rho(theta_t, tcorr_array);
        
        
        % Calculate/prepare outputs
        [q0, j0]        = obj.calcCharges(theta_0, c_idx(2), 0);
        CM              = zeros(obj.M, obj.M, 2, length(tcorr_array));
        theta           = cat(4, theta_0, theta_t);
        C1P             = zeros(obj.M, length(tcorr_array) + 1);
         
        for k = 1:length(tcorr_array)
            tic
            fprintf('Calculating correlation matrix (%d/%d) ... \n', k , length(tcorr_array));
            
            % Calculate g-functions, which will be acted on by the
            % correlation propagators.
            hi = obj.getOneParticleEV(c_idx(1), obj.rapid_grid); % should be (1xN)
            hj = obj.getOneParticleEV(c_idx(2), obj.rapid_grid);

            g_t = obj.applyDressing( hi, theta_t(:,:,:,k), tcorr_array(k) ); % should be (1xNxM)
            g_0 = obj.applyDressing( hj, theta_0, 0 );
            
            [qt, jt] = obj.calcCharges(theta_t(:,:,:,k), c_idx(1), tcorr_array(k));
            
            C1P(:,1)    = q0;
            C1P(:,1+k)  = qt;
            
            if areCurrents(1)
                g_t = g_t .* obj.calcEffectiveVelocities(theta_t(:,:,:,k), tcorr_array(k), obj.x_grid, obj.rapid_grid);
                C1P(:,1+k)  = jt;
            end
            if areCurrents(2)
                g_0 = g_0 .* obj.calcEffectiveVelocities(theta_0, 0, obj.x_grid, obj.rapid_grid);
                C1P(:,1)    = j0;
            end
            
            
            % Calculate correlation matrix
            CM(:,:,:,k) = GHDcorr.calcCorrelationMatrix(g_0, g_t, theta_t(:,:,:,k), rho_t(:,:,:,k), rhoS_t(:,:,:,k), u_t(:,:,:,k));
            
            toc
            fprintf('\n')
        end
    end
    
    
end % end public methods


methods (Access = protected)
    
    function Q_dr = applyDressing(obj, Q, theta, t)
        % Q must be defined on whole rapidity grid, but can depend on up to
        % two rapidities. i.e. dim(Q) = (GxNxM) or (GxNx1)
        % When performing dressing, it is most practical if second rapidity
        % argument is second dimension (might change this to be consistent
        % in all code in the future). Thus, within this method indices 1
        % and 2 are considered permuted!!!
        % NOTE: in correlation paper a distinction is made between vector
        % and scalar field. By these conventions T (diff scattering) is
        % symmetric, whereby scalar and vector fields are treated equally.
        
        % Check dimensionalities
        [d1, d2, d3] = size(Q);
        if d2 ~= obj.N 
            error('Quantity for dressing has wrong format!')
        end
        
        % Since we're working with permuted indices ...
        Q = permute(Q, [2 1 3]);
        
        % Set flags for homogeniety
        homo_state_flag     = 0; % is 1, if size(theta,3) == 1
        homo_quant_flag     = 0; % is 1, if size(Q,3) == 1
        if d3 == 1
            homo_quant_flag = 1;
        end
        if size(theta,3) == 1
            homo_state_flag = 1;
        end
        
        % Calculate dressing operator
        dk      = obj.rapid_grid(2) - obj.rapid_grid(1); 
        kernel  = dk/2/pi*obj.calcScatteringRapidDeriv(t, obj.x_grid, obj.rapid_grid, obj.rapid_grid); % is symmetric
        U       = eye(obj.N,obj.N) + kernel.*theta; % theta has been "transposed" - 2nd index is now 2nd rapid
        
        % We now have the equation Q = U*Q_dr. Therefore we solve for Q_dr
        % using the '\' operation.
        if homo_quant_flag && homo_state_flag % both are homogeneous
            Q_dr = U\Q;
        elseif homo_state_flag && ~homo_quant_flag % only Q is inhomogeneous
            Q_dr = zeros(d2, d1, d3);
            
            for i = 1:d3
                Q_dr(:,:,i) = U\Q(:,:,i);
            end
        elseif ~homo_state_flag && homo_quant_flag % only state is inhomogeneous
            Q_dr = zeros(d2, d1, size(theta,3));
            
            for i = 1:size(theta,3)
                Q_dr(:,:,i) = U(:,:,i)\Q;
            end
        else % both are inhomogeneous
            if d3 ~= size(theta,3)
                error('Quantity and theta have conflicting size!')
            end
            
            Q_dr = zeros(d2, d1, d3);
            
            for i = 1:d3
                Q_dr(:,:,i) = U(:,:,i)\Q(:,:,i);
            end
        end
        
        % Permute back to regular index convention
        Q_dr = permute(Q_dr, [2 1 3]);
           
    end
    
    
    function e_eff = calcEffectiveEnergy(obj, T, t, x, rapid)
        % Calculates the dressed energy per particle epsilon(k), from which the
        % preasure and later the filling factor theta can be derived.
        % This is achieved by iteratively solving the TBA equation.
        ebare       = obj.getBareEnergy(t, x, obj.rapid_grid); % shoud be (1xNxM) % MIGHT REPLACE BACK TO RAPID INPUT LATER
        dT          = obj.calcScatteringRapidDeriv(t, x, obj.rapid_grid, obj.rapid_grid );
        dk          = obj.rapid_grid(2) - obj.rapid_grid(1);
        
        % Need quantities in 3D format even if homogenous
        if size(ebare,3) == 1
            ebare = repmat(ebare, 1, 1, obj.M);
        end
        if size(dT,3) == 1
            dT = repmat(dT, 1, 1, obj.M);
        end
        
        e_eff       = zeros(1, obj.N, obj.M);
        e_eff_old   = zeros(1, obj.N, obj.M);
        error_rel   = 1;
        count       = 0;
        
        % Solve TBA eq. for epsilon(k) by iteration:
        % Using epsilonk_old, update epsilon_k until convergence is reached 
        % (difference between epsilonk_old and epsilon_k becomes less than tol)
        while any(error_rel > obj.tolerance) & count < obj.maxcount % might change this
            
            % calculate epsilon(k) from integral equation using epsilonk_old
            % i.e. update epsilon^[n] via epsilon^[n-1]
            for i = 1:obj.M
                % Transpose to column vector for matrix multiplication
                e_eff(:,:,i) = ebare(:,:,i)/T + transpose(dk/2/pi*dT(:,:,i)*log( 1 + exp(-transpose(e_eff_old(:,:,i))) ));
            end
            
            % calculate error
            sumeff      = sum( e_eff.^2 ,2);            
            error_rel   = squeeze(sum( (e_eff - e_eff_old).^2 )./sumeff);
            e_eff_old   = e_eff;

            count       = count+1;
        end
    end
    
    
    function [theta_next, theta_tmp, u_next] = performFirstOrderStep(obj, theta, theta_prev, u_prev, t, dt)            
        % Use single Euler step
        [v_eff, a_eff]  = obj.calcEffectiveVelocities(theta_prev, t, obj.x_grid, obj.rapid_grid); % should be (1xNxM)
            
        x_back          = obj.x_grid - dt*v_eff;
        r_back          = obj.rapid_grid - dt*a_eff;

        % Use interpolation to find theta_prev at x_back, r_back and
        % assign values to theta.
        theta_next      = obj.interpPhaseSpace(theta_prev, r_back, x_back, obj.extrapFlag);
        u_next          = obj.interpPhaseSpace(u_prev, r_back, x_back, true ); % always extrapolate u
        theta_tmp       = zeros(1,obj.N,obj.M); % Unused outout for similarity to higher order steps       
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
            [v_eff, a_eff] = obj.calcEffectiveVelocities(theta_mid, t+dt/2, obj.x_grid, obj.rapid_grid);

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
    
    
    function function_int = interpPhaseSpace(obj, function_grid, rapid_int, x_int, extrapFlag)
        % This function exists because MATLAB has different syntax between
        % interp1 and interp2, and I want to choose whether i extrapolate
        % or not with a simple TRUE/FAlSE argument.
        % ASSUME function_grid is on x_grid and rapid_grid.
        % Returns function_int with same dimensions as input function_grid
        
        % Need spacial dimension as first index in order to use (:) linearization
        x_int       = permute(x_int, [3 2 1]);
        rapid_int   = permute(rapid_int, [3 2 1]);
        
        % Determine whether function must be permuted for interpolation
        perm_flag = 0;
        if size(function_grid,3) == obj.M
            function_grid = permute(function_grid, [3 2 1]);
            perm_flag = 1;
        end
            
        
        if extrapFlag
            function_int = interp2( obj.rapid_grid, permute(obj.x_grid, [3 2 1]), function_grid, rapid_int(:), x_int(:), 'spline');
        else
            % Set all extrapolation values to zero!
            function_int = interp2( obj.rapid_grid, permute(obj.x_grid, [3 2 1]), function_grid, rapid_int(:), x_int(:), 'spline', 0);
        end
        
        % Returns output function to input dimensions
        function_int = reshape(function_int, obj.M, obj.N);
        if perm_flag
            function_int = permute(function_int, [3 2 1]);
        end
    end
  
end % end protected methods

end % end classdef