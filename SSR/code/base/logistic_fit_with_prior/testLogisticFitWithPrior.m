clear all
close all

% All vectors need to be column(!) vectors

ti = [1 2 3 6 7]';
yi = [0.1 0.15 0.5 0.75 0.8]';
yip = [0.4 0.475 0.5 0.525 0.55]';
c = 1;

gammas1 = [0.1, 0.1, 0.3, 0.1, 0.1]';  % this needs to be adapted (and modulated locally with 1/sigma(x)^2, one value per time-point!!
gammas2 = [10 10 10 10 10]';           % low variance (we trust the prior)
b0 = 0;
k0 = 0;

[bhat1,khat1] = logistic_fit_with_prior( ti, yi, yip, c, gammas1, b0, k0 );
[bhat2,khat2] = logistic_fit_with_prior( ti, yi, yip, c, gammas2, b0, k0 );

ls = linspace( 0, 10, 100 )';

lf1 = logistic_fcn( c, bhat1, khat1, ls );
lf2 = logistic_fcn( c, bhat2, khat2, ls );

% now plot it

figure
hold on
plot( ti, yi, 'b' )
plot( ti, yip, 'r' )
plot( ls, lf1, 'b-.' )
plot( ls, lf2, 'r-.' )

legend( 'data', 'prior data', 'fit (low gamma)', 'fit (high gamma)', 'location', 'best' );


