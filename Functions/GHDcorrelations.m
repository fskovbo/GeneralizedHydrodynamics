classdef GHDcorrelations
    
properties (Access = private)
    M           = [];
    N           = [];
    
    rapid_grid  = [];
    x_grid      = [];
       
    theta_0     = []; % filling function
    rho_0       = []; % density of particles
    rhoS_0      = []; % density of states
    a_eff0      = []; % effective acceleration at t = 0
    
    % Function handles copied from GHD
    applyDressing           = [];
    calcScatteringRapidDeriv= [];
    getStatFactor           = [];
    

end % end private properties


methods (Access = public)
    
    % constructor
    function obj = GHDcorrelations(x_grid, rapid_grid, theta_0, rho_0, rhoS_0, a_eff0, GHDfunctionhandles)
        % Check dimensionalites
        if ~isequal( size(theta_0), size(rho_0), size(rhoS_0), size(a_eff0) )
            error( 'State variables must have same dimensions!' )
        end
        if ~isequal( size(theta_0), [1, length(rapid_grid), length(x_grid)] )
            error( 'State variables do not match grids!' )
        end
        
        
        % Input GHDfunctionhandles as struct and copy fields into object properties.
        fn = fieldnames(GHDfunctionhandles);
        for l = 1:length(fn)                      
            if isprop(obj, fn{l}) % only copy field if defined among properties
                eval(['obj.',fn{l},' = GHDfunctionhandles.',fn{l},';']);
            end
        end
        
        % Copy input variables properties
        obj.N           = length(rapid_grid);
        obj.M           = length(x_grid);
        obj.x_grid      = x_grid;
        obj.rapid_grid  = rapid_grid;
        obj.a_eff0      = a_eff0;
        obj.theta_0     = theta_0;
        obj.rho_0       = rho_0;
        obj.rhoS_0      = rhoS_0;
    end
    

    function corrMat = calcCorrelationMatrix(obj, g_0, g_t, theta_t, rho_t, rhoS_t, u_t)
        % Maybe have time evolution here instead of in constructor?? would
        % kinda make more sense
        % Find better way to determine which type of correlation         
        
        corrMat     = zeros(obj.M, obj.M, 2);
        dudr_t      = gradient(u_t, obj.rapid_grid, 0, 0 ); % ***

        parfor y_index = 1:obj.M
%         for y_index = 1:obj.M

            theta_0y    = obj.theta_0(:,:,y_index);
            f_0y        = obj.getStatFactor(theta_0y);
            g_0y        = g_0(:,:,y_index);
            
            corr_prod    = rhoS_t.*theta_0y.*f_0y.*g_0y./abs(dudr_t); % (1xNxM)
         
            W2_temp    = obj.rhoS_0(:,:,y_index).*f_0y.*g_0y; % should be (1xNx1) 
            W2_temp_dr = obj.applyDressing( W2_temp, theta_0y, 0 ); % (1xNx1)
            W2_temp_sdr= W2_temp_dr - W2_temp; 
            
            W1         = 0;
            W3         = 0;
            
            corrMatCol      = zeros(obj.M,1,2);
            
            for x_index = 1:obj.M
                u_tx        = u_t(:,:,x_index);
                dudr_tx     = dudr_t(:,:,x_index);
                theta_tx    = theta_t(:,:,x_index);
                rho_tx      = rho_t(:,:,x_index);
                rhoS_tx     = rhoS_t(:,:,x_index);
                cp_tx       = corr_prod(:,:,x_index);
                g_tx        = g_t(:,:,x_index);
                
                gamma       = obj.findRootSet( u_tx, obj.x_grid(y_index), dudr_tx ); % ** (might give problems with extrapolation)

                
                % Calculate contribution from direct propagator
                direct_temp = obj.interpRapidities( cp_tx.*g_tx, gamma); % (1xG) ***
                direct      = sum(direct_temp,2); % sum over gamma ***
                
                % Calculate indirect propagator
                W2          = - heaviside( u_tx - obj.x_grid(y_index) ) .* W2_temp_sdr;
                [indirect, W1, W3] = obj.calcIndirectCorrelation(W1, W2, W3, gamma, theta_tx, rho_tx, rhoS_tx, u_tx, g_tx, cp_tx );
                
                corrMatCol(x_index,:,1) = direct;
                corrMatCol(x_index,:,2) = indirect;
                
            end
            corrMat(:,y_index,:) = corrMatCol;
        end
        
    end
    
end % end public methods


methods (Access = private)

    function [indirect, W1, W3] = calcIndirectCorrelation(obj, W1, W2, W3, gamma, theta_tx, rho_tx, rhoS_tx, u_tx, g_tx, corr_prod_x )
        % Moved calculation to separate function to declutter for-loops
        
        if all( obj.a_eff0 == 0 ) % homogeneous system
            indirect    = 0;
            return
        end
        
        dk      = obj.rapid_grid(2) - obj.rapid_grid(1);
        dx      = obj.x_grid(2) - obj.x_grid(1);
        f_tx    = obj.getStatFactor( theta_tx );

        
        % Update W1 
        Kern    = 1/(2*pi) * obj.calcScatteringRapidDeriv( obj.rapid_grid, gamma );
        Kern_dr = obj.applyDressing(Kern, theta_tx, 1);
        W1      = W1 - dx*Kern_dr.*obj.interpRapidities(corr_prod_x , gamma); 

        
        % Solve for Delta
        a_eff   = obj.interpSpace( obj.a_eff0, u_tx ); % evaluate a_eff0 at x = u(t_corr, x_corr, lambda)         
        a_eff   = permute( diag( permute(a_eff, [3 2 1])), [2 1]); % a_eff and u must be evaluated at same rapidity  

        kernel  = dk/(2*pi) * obj.calcScatteringRapidDeriv( obj.rapid_grid, obj.rapid_grid );
        U       = eye(obj.N,obj.N) + kernel.*theta_tx;

        vec     = rhoS_tx .* f_tx;
        Umat    = diag( (2*pi*a_eff).^(-1) ) + dx*diag( vec ) - dx*inv(U).*vec; % is right, if U is right --- otherwise use U' 

        Delta   = (Umat\(W1+W2+W3)')'; % Solve integral equation through iteration


        % IMPORTANT! update W3 for next step
        integr      = vec .* Delta;
        integr_sdr  = obj.applyDressing( integr, theta_tx, 1) - integr; % should be (1xNxM)
        W3          = W3 + dx*integr_sdr; 

        
        % Calculate indirect contribution via Delta
        indirect    = trapz( obj.rapid_grid, Delta.*rho_tx.*f_tx.*g_tx ,2); % integrate over rapidity
    end
    
    
    function gamma = findRootSet(obj, u_xt, y, du_xt)
        % Finds gamma such that u(x,t,gamma) = u_xt(gamma) = y
        
        if nargin < 4
            du_xt = gradient(u_xt, obj.rapid_grid, 0);
        end
        
        % Turn u_xt into continuous, anonymous function for finding roots
        % NOTE: If fzero fails, it is often due to extrapolation. Thus, Try
        % different algorithms for interp1!
%         u_func = @(rapid) pchip(obj.rapid_grid, u_xt, rapid);
        u_func = @(rapid) interp1(obj.rapid_grid, u_xt, rapid, 'linear','extrap');
        
        if all( du_xt < 0)
            % If du_xt is monotonically decreasing with rapidity, the root
            % set contains only a single member.
            
            gamma = fzero(@(x) u_func(x) - y, 0);
        else
            % Multiple members in root set. Use sign flip to gauge how
            % many. 
            % NOTE: does not take into account u = y continuously
            
            signflip    = diff( u_xt - y >=0 ) ~= 0; % logical N-1 length vector indicating signflips
            rapid0      = obj.rapid_grid( [signflip false] );
            Nroots      = length(rapid0);
            gamma       = zeros(1,Nroots);
            
            for i = 1:Nroots
                gamma(i) = fzero(@(x) u_func(x) - y, rapid0(i));
            end    
            
            gamma       = unique(gamma);
        end
    end
    
    
    function out = interpRapidities(obj, table, lookupvals)
        % interp1 interpolates each column if ndims > 1
        % make sure rapidity is first index
        table = permute(table, [2 3 1]);

        out = interp1(obj.rapid_grid, table, lookupvals, 'spline');
        % out has dimensions (GxM), where G is length(lookupvals) and M is
        % spatial dimension (if included in input).
        % We want typical outpu format of (1xGxM)

        out = permute(out, [3 1 2]);
    end
    
    
    function out = interpSpace(obj, table, lookupvals)
 
        if ndims(table) < 3
             error('Quantity for spatial interpolation is wrong format!')
        end
        
        % interp1 interpolates each column if ndims > 1
        % make sure space is first index
        table = permute(table, [3 2 1]);

        out = interp1(permute(obj.x_grid, [3 2 1]), table, lookupvals, 'spline');
        % out has dimensions (NxK), where K is length(lookupvals) and N is
        % rapidity dimension.
        % We want typical output format of (1xNxK)

        out = permute(out, [3 2 1]);
    end
    
end % end private methods
    
end % end classdef 
