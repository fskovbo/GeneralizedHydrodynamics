function rho_fit = calcZeroTempState(x_grid_si, k_grid_si, c_si, U_si, Natoms)
    
    % Convert to internal units
    m_si    = 87*1.6605402e-27; % Rb-87 mass
    hbar_si = 1.054571726e-34;

    Eg_si   = 0.5 * hbar_si^2 * c_si^2 / m_si; 
    Lg_si   = 1/c_si;
    
    x_grid  = x_grid_si/Lg_si;
    k_grid  = k_grid_si*Lg_si;
    U       = U_si/Eg_si;
    
    
    % Fit mu0 to Natoms
    mu0_guess   = U(1);
    fitfunc     = @(mu0) calcAtomnumberDiff(mu0, x_grid, k_grid, U, Natoms);
    options     = optimset('Display','iter');
    mu0_fit     = fminsearch(fitfunc, mu0_guess, options);
    
    rho_fit     = calcRho( mu0_fit, k_grid, U );
end


function f = solveF(k_grid)
    N       = length(k_grid);
    dk      = k_grid(2) - k_grid(1);

    f       = zeros(1,N);
    f_old   = zeros(1,N);

    error_rel   = 1;
    count       = 0;
    maxcount    = 1e4;
    tolerance   = 1e-8;

    while (error_rel > tolerance & count < maxcount) 

        % Calculate g(y)
        f           = 1/(2*pi) + dk/pi*sum( f_old./(1 + (k_grid-k_grid').^2) , 2 );

        % Calculate error
        sumg        = sum( f.^2 );            
        error_rel   = sum( (f - f_old).^2 )./sumg;
        f_old       = f;

        count       = count+1;
    end

end


function rho = calcRho( mu0, k_grid, U )
    mu      = max( mu0 - U, 1e-8); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rho     = zeros( length(k_grid), length(mu) );    
    
    for i = 1:length(mu)
        K       = sqrt(mu(i));
        K_grid  = linspace(-K, K, length(k_grid));
        f       = solveF(K_grid);
        
        rho(:,i)= interp1(K_grid, f, k_grid);
    end
    
    rho( isnan(rho) ) = 0; % !!!!!!!!!!!!!!!!!!!!!!
end


function Ndiff = calcAtomnumberDiff(mu0, x_grid, k_grid, U, Natoms)
    rho     = calcRho( mu0, k_grid, U );    
    Ncalc   = trapz(k_grid, trapz(x_grid, rho'));
    
    Ndiff   = abs( Natoms - Ncalc );
end
