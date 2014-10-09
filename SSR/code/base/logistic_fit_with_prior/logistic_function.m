function ret = logistic_function(x,a,b,varargin)
% Computes the logistic function
%
% F(x;P) = 1/(1+exp(-1/B*(x-A))
%

A = 0;
K = 1;

% use glmfit parameters
if( numel(varargin)>0 ), A = varargin{1}; end
if( numel(varargin)>1 ), K = varargin{2}; end

if( numel(varargin)>2 )
     b = 1/b;
     a = -a*b;
end

%ret = 1./(1+exp(-1/b*(x-a)));

ret = A + (K-A)./(1+exp(-1/b*(x-a)));

%% glmfit vs logistic_function
% xx = linspace(-5,5,100);
% 
% b = @(x0,x1) (x1+x0)/2;
% a = @(x0,x1,y) -log(1/y - 1)/(x1-b(x0,x1));
% y = 0.999;
% 
% x0 = -0.5; x1 = 3;
% yy1 = glmval([-b(x0,x1)*a(x0,x1,y) a(x0,x1,y)]', xx, 'logit' );
% yy2 = logistic_function(xx,b(x0,x1),1/a(x0,x1,y));
% plot(xx,yy1,'b-',xx,yy2,'r--',x0,1-y,'go',x1,y,'go'); grid on;
