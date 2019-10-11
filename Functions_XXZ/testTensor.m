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

%%
% 
% lala = GHDtensor( reshape(1:24, 3, 4, 1, 1, 2) );
% 
% yoo2 = atSpacePoint(lala,1)

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