clear all; clc;

P = '2*atan( coth(t*T/2) * tan(k) )';
D = '1.5+0.3*tanh(3*t).*sin(4*pi*(x-t))';

dP = takeDerivStr(P, 'k')

dD = takeDerivStr(D,'x') 


function [deriv_str, isZero] = takeDerivStr(func_str, var_str)
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
    
    deriv_sym   = simplify( deriv_sym , 'Steps', 100);

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