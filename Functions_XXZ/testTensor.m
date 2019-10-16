clear all; clc;

hej = reshape(1:24, 3, 4 , 2)

test = GHDtensor(hej);

yoo1 = test(:,:,1)

test(1,1,1) = 5;

t{1} = test.*hej;
t{2} = test.*test;
t{3} = hej.*test;
t{4} = test*0.5;
t{5} = 0.5*test;
t{6} = -test;
t{7} = test + test;
t{8} = test + 1;
t{9} = 1+test;
t{10} = test - test;
t{11} = test - 1;
t{12} = 1 - test;
t{13} = exp(test);
t{14} = log(test);
t{15} = log(exp(test));
t{16} = test./4;
t{17} = 4./test;
t{18} = test/4;
t{19} = 4/test;
t{20} = test./hej;


for i = 1:length(t)
    disp(' ========================================== ')
    double(t{i})
end

%% Test Multiplication

A = sym('A', [4 3 2 2]);
B = sym('B', [3,1,2]);

At = GHDtensor(A);
Bt = GHDtensor(B);

x = At*Bt;

% in matlab sum only works in 2D for symbolic
% x_lm^in = sum_jk A_lk^ij * B_km^jn

for l = 1:size(A,1)
for m = 1:size(B,2)
for i = 1:size(A,3)
for n = 1:size(B,4)
    for j = 1:size(A,4)
        temp(j) =  A(l,:,i,j) * B(:,m,j,n); % sum over k
    end
    x_ref(l,m,i,n) = sum( temp ); % sum over j
end
end
end
end

isequal(x_ref, double(x))


%% Test Linear solve
% A = sym('A', [4 3 2 2]);
% B = sym('B', [4 2 2 1]);
A = rand(4,5,1,1 );
B = rand(4,2,1,1 );

At = GHDtensor(A);
Bt = GHDtensor(B);

x = At\Bt;

% in matlab sum only works in 2D for symbolic
% Check that A*x = B

compare = ismembertol( double(At*x) , double(Bt), 1e-10 );
all( compare(:) )


%% Test inv()

A   = rand(4,4,2,2 );
B   = rand(4,2,2,1 );

At  = GHDtensor(A);
Bt = GHDtensor(B);

Ai  = inv(At);

double(At*Ai)

x1  = double(At\Bt)
x2  = double(Ai*Bt)

%%
% 
% U = reshape( 1:3*3*2*2, 3, 3, 2, 2 );
% 
% U = rand(3,3,2,2,2);
% Ut = GHDtensor(U)
% 
% % h = reshape( 1:3*2, 3, 1, 2, 1 );
% h = rand(3,1,2,1);
% ht = GHDtensor(h)
% 
% 
% double(Ut\ht)