function tform = tform_resize( tform, m )

if( ~isequal(tform.m, m) )
    ndim = numel( m );
    yc = reshape( tform.yc, [tform.m ndim] );
    xc = getCenteredGrid( tform.omega, m );

    yct = zeros( [m ndim] );
    
    for i=1:ndim
        switch( ndim )
            case 2, yct(:,:,i) = reshape( inter(yc(:,:,i), tform.omega, xc), m );                
            case 3, yct(:,:,:,i) = reshape( inter(yc(:,:,:,i), tform.omega, xc), m  );
        end
        
    end
    
    tform.yc = reshape( yct, [], 1 );
    tform.m  = m;    
end