classdef testFormat
   
properties (Access = public)
    couplingNames = {'mu' , 'c'}
    energy = 'rapid^2 - mu'
    momentum = 'rapid'
    scattering = '-2*atan(rapid/c)'
    
    couplings = []
    mainDerivs = []
    couplingDerivs = []
end


methods (Access = public)
    
    function obj = testFormat(couplings)
        
        obj.couplings = couplings;
        
        dedr = obj.takeDerivStr( obj.energy, 'rapid' );
        dpdr = obj.takeDerivStr( obj.momentum, 'rapid' );
        dTdr = obj.takeDerivStr( obj.scattering, 'rapid' );
        
        obj.mainDerivs{1,1} = obj.formatFunction(dedr, false);
        obj.mainDerivs{2,1} = obj.formatFunction(dpdr, false);
        obj.mainDerivs{3,1} = obj.formatFunction(dTdr, true);
        
        for i = 1:length(obj.couplingNames)
            dedc = obj.takeDerivStr( obj.energy, obj.couplingNames{i} );
            dpdc = obj.takeDerivStr( obj.momentum, obj.couplingNames{i} );
            dTdc = obj.takeDerivStr( obj.scattering, obj.couplingNames{i} );

            obj.mainDerivs{1,i+1} = obj.formatFunction(dedc, false);
            obj.mainDerivs{2,i+1} = obj.formatFunction(dpdc, false);
            obj.mainDerivs{3,i+1} = obj.formatFunction(dTdc, true);
        end
        
        
        for i = 1:length(obj.couplingNames)           
            coup_str = func2str( obj.couplings{i} );
            coup_str = erase(coup_str, ' '); % erase spaces
            coup_str = erase(coup_str, '@(t,x)'); % erase function prefix
            dcdx = obj.takeDerivStr( coup_str, 'x' );
            dcdt = obj.takeDerivStr( coup_str, 't' );
            
            obj.couplingDerivs{1,i} = obj.formatFunction(dcdt, false);
            obj.couplingDerivs{2,i} = obj.formatFunction(dcdx, false);
        end
    end
    
    
    function deriv_str = takeDerivStr(obj, func_str, var_str)
        % Takes a function formated as a string and output its derivative
        % (with respect to the string variable var_str) as a string
        
        % Convert to symbolics
        func_sym    = str2sym( func_str );
        var_sym     = str2sym( var_str );
        deriv_sym   = diff( func_sym, var_sym );
        
        % If derivative is zero, return empty 
        if isequal(deriv_sym, sym(0))
            deriv_sym = [];
        end
        
        % Format back to string
        deriv_str   = char( deriv_sym );
    end
    
    
    function func = formatFunction(obj, func_str, isKernel)
        % Takes a function formated as a string and outputs an anonymous function 
        
        if isempty(func_str)
            if isKernel
                % kernel should not be empty
                func_str = '0'; 
            else
                func = [];
                return
            end
        end
        
        % Replace operations with elementwise ones
        func_str = regexprep(func_str, '*', '.*');
        func_str = regexprep(func_str, '/', './');
        func_str = regexprep(func_str, '\^', '\.^');
        
        % Replace all couplings with their expressions
        for i = 1:length(obj.couplingNames)
            coup_str = func2str( obj.couplings{i} );
            coup_str = erase(coup_str, ' '); % erase spaces
            coup_str = erase(coup_str, '@(t,x)'); % erase function prefix
            func_str = regexprep(func_str, obj.couplingNames{i}, coup_str);
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
        func     = str2func(func_str);
    end
end


end