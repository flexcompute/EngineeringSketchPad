function adjoint2
    % written by John Dannenhoffer

    % design parameters
    syms xbeg ybeg zbeg xlen ylen zlen radius
    
    % mass properties
    area = 2 * xlen * ylen + 2 * ylen * zlen + 2 * zlen * xlen ...
         - 2 * pi * radius^2 + 4 * pi * radius * zlen
     
    vol = xlen * ylen * zlen
     
    xcg = ((xbeg+xlen/2) * xlen * ylen * zlen + (xlen /   3) * pi * radius^2 * zlen) / vol
     
    ycg = ((ybeg+ylen/2) * xlen * ylen * zlen + (ylen   / 3) * pi * radius^2 * zlen) / vol
     
    zcg = ((zbeg+zlen/2) * xlen * ylen * zlen + (zlen      ) * pi * radius^2 * zlen) / vol
    
    % loop through the mass properties
    for foo = [vol, area, xcg, ycg, zcg]
    
        % print nominal value
        val = double(subs(foo, {xbeg, ybeg, zbeg, xlen, ylen, zlen, radius}, ...
                               {-1,   -2,    -3,  4,    3,    1,    2/3    }));
        fprintf(1, '%s = %12.4f\n', foo, val);
    
        % print derivatives wrt the design parameters
        for bar = [xbeg, ybeg, zbeg, xlen, ylen, zlen, radius]
            val = double(subs(diff(foo,bar), {xbeg, ybeg, zbeg, xlen, ylen, zlen, radius}, ...
                                             {-1,   -2,    -3,  4,    3,    1,    2/3    }));
                                   
            fprintf(1, "d(%s)/d(%s) = %12.4f\n", foo, bar, val);
        end % for bar
    end % for foo
           
end % function adjoint1
     
     
     
