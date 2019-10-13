clearvars
clc

p = sym('p', [20,1]);
p(99) = 0;

Ntypes = 5;
type1 = 1:Ntypes;
type2 = 1:Ntypes;

type1   = reshape(type1, Ntypes,1); % type1 is 1st index
type2   = reshape(type2, 1, Ntypes); % type2 is 2nd index


% Calculate 3 term by creating meshgrid type index matrices 
m1      = repmat(type1, 1, Ntypes);
m2      = repmat(type2, Ntypes, 1);
min_idx = min(m1,m2) - 1; % matrix of minimum index

dT3 = 0;
D = abs(type1-type2);
type_arg = zeros(Ntypes);
for n = 1:Ntypes-1
    % dT3^{i,j} ~ sum (n to min(i,j)-1) dp^{ |i-j| + 2*n }
    logi     =  min_idx >= n;
    type_arg = D.*logi + 2*n*logi;
    dT3      = dT3 + getP(p,type_arg);
end

dT4 = sym(zeros(Ntypes,Ntypes));
for i = 1:Ntypes
    for j = 1:Ntypes
        for n = (abs(i-j)+2):2:(i+j-2)
            dT4(i,j) = dT4(i,j) + getP(p,n);
        end
    end
end



dT3
dT4


%%
clear all; clc
syms rapid Delta type 

p = 2*atan( coth( type.*acosh(Delta)/2 ) .* tan(rapid) )

% diff(p, Delta)
deriv = simplify(diff(p, rapid))

deriv_str = char(deriv)
deriv_str = regexprep(deriv_str,'*','.*')
fun_str = ['@(rapid, type) ' deriv_str] 

my_fun = str2func(fun_str)

%%

couplings{1} = @(t,x) 5 - x.^2;
couplings{2} = @(t,x) t;


test = testFormat(couplings)

test.mainDerivs
test.couplingDerivs



%%
% 
% function out = getP(p, index)
%     index( index < 1 ) = 99;
%     out = p(index);
% end
