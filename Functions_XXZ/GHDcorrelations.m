classdef GHDcorrelations
    % Class containing necessary functions for calculating GHD correlation
    % functions. Has been separated from other classes to avoid clutter.
    % NOTE: Currently only supports initial states (marked with _0) which
    % are TBA stationary states. Furthermore, time evolved states (marked
    % with _t) must be homogeneous propagations of initial state!
    
    
properties (Access = private)
    M           = [];
    N           = [];
    Ntypes      = [];
    
    rapid_grid  = [];
    x_grid      = [];
       
    theta_0     = []; % filling function
    rho_0       = []; % density of particles
    rhoS_0      = []; % density of states
    a_eff0      = []; % effective acceleration at t = 0
    
    % Function handles copied from GHD solver
    applyDressing           = [];
    calcScatteringRapidDeriv= [];
    getStatFactor           = [];
    

end % end private properties


methods (Access = public)
    
    % constructor
    function obj = GHDcorrelations(x_grid, rapid_grid, theta_0, rho_0, rhoS_0, a_eff0, GHDfunctionhandles)       
        
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
        obj.Ntypes      = size(theta_0,3);
        obj.x_grid      = x_grid;
        obj.rapid_grid  = rapid_grid;
        obj.a_eff0      = a_eff0;
        obj.theta_0     = theta_0;
        obj.rho_0       = rho_0;
        obj.rhoS_0      = rhoS_0;
    end
    

    function corrMat = calcCorrelationMatrix(obj, g_0, g_t, theta_t, rho_t, rhoS_t, u_t)
        % Calculates correlations between g_0 and g_t
        
        corrMat     = zeros(obj.M, obj.M, 2);
        [~, dudr_t] = gradient(double(u_t), 0, obj.rapid_grid, 0, 0, 0 ); %
        dudr_t      = GHDtensor( dudr_t );

%         parfor y_index = 1:obj.M
        for y_index = 1:obj.M

            theta_0y    = atSpacePoint(obj.theta_0, y_index);
            f_0y        = obj.getStatFactor(theta_0y);
            g_0y        = atSpacePoint(g_0, y_index);
            
            corr_prod    = rhoS_t.*theta_0y.*f_0y.*g_0y./abs(dudr_t); 
         
            % Calculate source terms required for indirect propagator
            W2_temp    = atSpacePoint(obj.rhoS_0,y_index).*f_0y.*g_0y;  
            W2_temp_dr = obj.applyDressing( W2_temp, theta_0y, 0 ); 
            W2_temp_sdr= W2_temp_dr - W2_temp; 
            
            W1         = 0;
            W3         = 0;
            
            corrMatCol      = zeros(obj.M,1,2);
            
            for x_index = 1:obj.M
                u_tx        = atSpacePoint(u_t, x_index);
                dudr_tx     = atSpacePoint(dudr_t, x_index);
                theta_tx    = atSpacePoint(theta_t, x_index);
                rho_tx      = atSpacePoint(rho_t, x_index);
                rhoS_tx     = atSpacePoint(rhoS_t, x_index);
                cp_tx       = atSpacePoint(corr_prod, x_index);
                g_tx        = atSpacePoint(g_t, x_index);
                
                gamma       = obj.findRootSet( u_tx, obj.x_grid(y_index), dudr_tx );

                
                % Calculate contribution from direct propagator
                direct_temp = obj.interp2Gamma( cp_tx.*g_tx, gamma); % (1xG)
                direct      = sum(double(direct_temp) ,1); % sum over gamma 
                
                % Calculate indirect propagator
                W2          = - heaviside( double(u_tx) - obj.x_grid(y_index) ) .* W2_temp_sdr;
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
        
        if all( double(obj.a_eff0) == 0 ) % homogeneous system --> no indirect corr
            indirect    = 0;
            return
        end
        
        delta_k = obj.rapid_grid(2) - obj.rapid_grid(1);
        delta_x = obj.x_grid(2) - obj.x_grid(1);
        f_tx    = obj.getStatFactor( theta_tx );

        
        % Update First source term, W1 
        Kern    = 1/(2*pi) * obj.calcScatteringRapidDeriv( obj.rapid_grid, gamma );
        Kern_dr = obj.applyDressing(Kern, theta_tx, 1);
        W1      = W1 - delta_x*Kern_dr.*obj.interp2Gamma(corr_prod_x , gamma); 

        
        % Solve for Delta
        a_eff = obj.interp2u( obj.a_eff0, u_tx ); % evaluate a_eff0 at x = u(t_corr, x_corr, lambda)
        
        kernel  = delta_k/(2*pi) * obj.calcScatteringRapidDeriv( obj.rapid_grid, obj.rapid_grid );
        I_rapid = eye(obj.N);
        I_type  = repmat(eye(obj.Ntypes), 1 ,1, 1, 1);
        I_type  = permute(I_type, [3 4 1 2]);
        identity= I_rapid.*I_type;

        U       = identity + kernel.*transpose(theta_tx); 
        Uinv    = inv(U);
        
        vec     = rhoS_tx .* f_tx;
        Tmat    = identity.*((2*pi*a_eff).^(-1) + delta_x*vec) - delta_x*Uinv.*transpose(vec); % vec should be transposed for consistency!!

        Delta   = Tmat\(W1+W2+W3); % Solve integral equation through iteration


        % IMPORTANT! update W3 for next step
        integr      = vec .* Delta;
        integr_sdr  = obj.applyDressing( integr, theta_tx, 1) - integr; % should be (1xNxM)
        W3          = W3 + delta_x*integr_sdr; 

        
        % Calculate indirect contribution via Delta
        indirect    = sum(trapz( obj.rapid_grid, double(Delta.*rho_tx.*f_tx.*g_tx) ,1), 3); % integrate over rapidity and sum over type
    end
    
    
    function gamma = findRootSet(obj, u_xt, y, du_xt)
        % Finds gamma such that u(x,t,gamma) = u_xt(gamma) = y
        % Note: gamma is matrix of size (G,1,Nt)
        
        if nargin < 4
            [~, du_xt] = gradient(double(u_xt), 0 , obj.rapid_grid);
        end
        
        u_xt    = double(u_xt);
        du_xt   = double(du_xt);
        gamma   = [];
        
        for i = 1:obj.Ntypes
            % Turn u_xt into continuous, anonymous function for finding roots
            % NOTE: If fzero fails, it is often due to extrapolation. Thus, Try
            % different algorithms for interp1!
    %         u_func = @(rapid) pchip(obj.rapid_grid, u_xt, rapid);
            u_func = @(rapid) interp1(obj.rapid_grid, u_xt(:,:,i), rapid, 'linear','extrap');

            if all( du_xt(:,:,i) < 0)
                % If du_xt is monotonically decreasing with rapidity, the root
                % set contains only a single member.

                gamma_i = fzero(@(x) u_func(x) - y, 0);
            else
                % Multiple members in root set. Use sign flip to gauge how
                % many. 
                % NOTE: does not take into account u = y continuously

                signflip    = diff( u_xt(:,:,i) - y >=0 ) ~= 0; % logical N-1 length vector indicating signflips
                rapid0      = obj.rapid_grid( [signflip false] );
                Nroots      = length(rapid0);
                gamma_i       = zeros(1,Nroots);

                for i = 1:Nroots
                    gamma_i(i) = fzero(@(x) u_func(x) - y, rapid0(i));
                end    

                gamma_i       = unique(gamma_i);
            end
            
            gamma(:,1,i) = gamma_i;
        end
    end
    
    
    function tensor_out = interp2Gamma(obj, tensor_in, gamma)
        % gamma should be matrix of dimesions (G,1,Nt,1,M) 
        % NOTE: interp1 interpolates each column if ndims > 1
        
        out_size    = size(tensor_in);
        out_size(1) = size(gamma,1);
        mat_out     = zeros(out_size);
        
        for i = 1:obj.Ntypes
            mat_out(:,1,i,1,:) = interp1(obj.rapid_grid, tensor_in(:,:,i,:,:), gamma(:,:,i), 'spline');
        end

        if numel(mat_out) == 1 % mat_out is scalar
            tensor_out      = GHDtensor(1,1,1,1,1);
            tensor_out(1)   = mat_out;
        else
            tensor_out      = GHDtensor(mat_out);
        end
    end
    
    
    function tensor_out = interp2u(obj, tensor_in, u_x)
        % Interpolates tensor spatially to T( u(x, lambda), lambda ),
        % where u_x(lambda) = u(x, lambda)
        % NOTE: interp1 interpolates each column if ndims > 1

        mat_out = zeros(obj.N, 1, obj.Ntypes, 1, 1);
        
        for i = 1:obj.Ntypes
            tens_i  = tensor_in(:,:,i,:,:);
            ux_i    = u_x(:,:,i);
            
            % make sure spatial index is first
            tens_i  = permute(tens_i, [ 5 2 3 4 1]);
            x_g     = permute(obj.x_grid, [5 2 3 4 1]);
            
            mat_int = interp1(x_g, tens_i, ux_i, 'spline'); % (N,1,1,1,N)
            
            % Since u(x, lambda) must be evaluated at same rapidities as
            % the tensor T, one has to take the diagonal.
            mat_int = diag( permute(mat_int, [ 1 5 3 4 2]) ); % (N,1,1,1,1)
            
            mat_out(:,:,i,:,:) = mat_int;
        end

        tensor_out = GHDtensor(mat_out);
    end
    
end % end private methods
    
end % end classdef 
