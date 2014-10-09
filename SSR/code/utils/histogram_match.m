function [im,stretch] = histogram_match( im1, im2, mask1, mask2, mask, varargin )
% Histogram normalization 
%
% im1                   : reference image
% im2                   : image to normalize
% mask1                 : mask used to compute im1 histogram
% mask2                 : mask used to compute im2 histogram
% mask                  : mask for applying normalization in im1
%
% OPTIONS:
%  'show_plot',<0|1>    : 0=no plotting (default), 1=plot histograms

show_plot = 0;

for k=1:2:length( varargin )
    eval( [ varargin{k}, '=varargin{', int2str(k+1), '};' ] );
end

% The masks have to be equal size -> sample the larger mask with
% N number of points, where N = size of smaller mask.

idx1 = find( mask1>0 ); n1 = numel( idx1 );
idx2 = find( mask2>0 ); n2 = numel( idx2 );

if( n1 ~= n2 )
    diff = abs( n1-n2 );
    % since randn number are not necessarily unique we might zero
    % out certain pixels multiple times, to avoid this:        
    ridx = randperm( max(n1,n2) );
    ridx = ridx(1:diff);
    if( n1 > n2 )
        mask1( idx1(ridx) ) = 0;
    else
        mask2( idx2(ridx) ) = 0;
    end                
end

% Normalize im2 to match im1
im2 = im2 - min(im2(:));
im2 = im2 * ( (max(im1(:))-min(im1(:))) / max(im2(:)) )  -  min(im1(:));


[hh,xx] = hist([im1(mask1),im2(mask2)],128); 
[~,hh_max] = max(hh);

% get mean and mode from the joint histogram
hh1_mean = sum( hh(:,1)/sum(hh(:,1)).*xx ); % mean = expected value
hh1_mode = xx(hh_max(1));

% img
hh2_mean = sum( hh(:,2)/sum(hh(:,2)).*xx ); 
hh2_mode = xx(hh_max(2));

stretch = [hh1_mode/hh2_mode, hh1_mean-hh1_mode*hh2_mean/hh2_mode];

im = im2;
im(mask) = stretch(1) * im(mask) + stretch(2);

if( show_plot )
    hhs = hist( im(mask2), xx );
    plot( xx, hh(:,1), xx, hh(:,2), xx*stretch(1)+stretch(2), hh(:,2), xx, hhs ); 
    legend('image 1','image 2','streched hist 2', 'image 2 corr')
end