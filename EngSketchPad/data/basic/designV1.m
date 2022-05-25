function designV1(ivar, iplot)
    % designV1 - demonstration of geometric and tessellation sensitivities
    %
    % inputs:
    %    ivar  =1 for len
    %          =2 for dip
    %          =3 for ang
    %          =4 for xA
    %          =5 for yA
    %    iplot =0 plot outline with corner sensitivities
    %          =1 superimpose geometric    sensitivities
    %          =2 superimpose tesselation  sensitivities
    %
    % written by John Dannenhoffer

    if (nargin < 2)
        iplot = 99;
    end % if

    % inputs
    len = 2;
    dip = 0.5;
    ang = 40 * (pi/180);
    xA  = 0.5;
    yA  = -.5;

    len_dot = 0;
    dip_dot = 0;
    ang_dot = 0;
    xA_dot  = 0;
    yA_dot  = 0;

    % user-specified sensitivity
    if     (ivar == 1)
        despmtr = 'len';
        len_dot = 1;
    elseif (ivar == 2)
        despmtr = 'dip';
        dip_dot = 1;
    elseif (ivar == 3)
        despmtr = 'ang';
        ang_dot = 1;
    elseif (ivar == 4)
        despmtr = 'xA';
        xA_dot  = 1;
    elseif (ivar == 5)
        
        despmtr = 'yA';
        yA_dot  = 1;
    end % if

    % overall (nodal) geometry
    xB     = xA     + len;
    xB_dot = xA_dot + len_dot;

    yB     = yA;
    yB_dot = yA_dot;

    xC     = xA;
    xC_dot = xA_dot;

    yC     = yA     + len     * tan(ang);
    yC_dot = yA_dot + len_dot * tan(ang) + ang_dot * len / cos(ang)^2;

    % arc AB
    R     = ((len/2)^2 + (dip)^2) / (2*dip);
    R_dot = ((2*dip    ) * (len/2*len_dot + 2*dip*dip_dot) ...
            -(2*dip_dot) * (len/4*len     +   dip*dip    )) / (2*dip)^2;

    xD     = xA     + len     / 2;
    xD_dot = xA_dot + len_dot / 2;

    yD     = yA     + R     - dip;
    yD_dot = yA_dot + R_dot - dip_dot;

    t1beg  = atan2(yA-yD, xA-xD);
    t1end  = atan2(yB-yD, xB-xD);

    tt     = linspace(t1beg, t1end, 51);

    x1     = xD     + R     * cos(tt);
    x1_dot = xD_dot + R_dot * cos(tt);

    y1     = yD     + R     * sin(tt);
    y1_dot = yD_dot + R_dot * sin(tt);

    nx1    = cos(tt);
    ny1    = sin(tt);

    g1     = x1_dot .* nx1 + y1_dot .* ny1;

    % line BC
    t2beg  = 0;
    t2end  = sqrt((xC-xB)^2 + (yC-yB)^2);

    tt     = linspace(t2beg, t2end, 51);
    x2     = xB     + (xC     - xB    ) * (tt - t2beg) / (t2end - t2beg);
    x2_dot = xB_dot + (xC_dot - xB_dot) * (tt - t2beg) / (t2end - t2beg);

    y2     = yB     + (yC     - yB    ) * (tt - t2beg) / (t2end - t2beg);
    y2_dot = yB_dot + (yC_dot - yB_dot) * (tt - t2beg) / (t2end - t2beg);

    nx2    = (yC - yB) / (t2end - t2beg);
    ny2    = (xB - xC) / (t2end - t2beg);

    g2     = x2_dot .* nx2 + y2_dot .* ny2;

    % line CA
    t3beg  = 0;
    t3end  = sqrt((xA-xC)^2 + (yA-yC)^2);

    tt     = linspace(t3beg, t3end, 51);

    x3     = xC     + (xA     - xC    ) * (tt - t3beg) / (t3end - t3beg);
    x3_dot = xC_dot + (xA_dot - xC_dot) * (tt - t3beg) / (t3end - t3beg);

    y3     = yC     + (yA     - yC    ) * (tt - t3beg) / (t3end - t3beg);
    y3_dot = yC_dot + (yA_dot - yC_dot) * (tt - t3beg) / (t3end - t3beg);

    nx3    = (yA - yC) / (t3end - t3beg);
    ny3    = (xC - xA) / (t3end - t3beg);

    g3     = x3_dot .* nx3 + y3_dot .* ny3;

    % start a plot with the 3 edges
    plot(x1, y1, 'b-', x2, y2, 'b-', x3, y3, 'b-')
    xlabel('x')
    ylabel('y')
    title(sprintf('sensitivity with respect to \"%s\"', despmtr))
    axis equal

    % superimpose "correct" solution at the nodes
    hold on
        plot([xA+xA_dot, xB+xB_dot, xC+xC_dot], ...
             [yA+yA_dot, yB+yB_dot, yC+yC_dot], 'ro')
    hold off

    % compute and show the geometric sensitivities
    if (iplot >= 1)
        x1geom = x1 + g1 .* nx1;
        y1geom = y1 + g1 .* ny1;

        x2geom = x2 + g2 .* nx2;
        y2geom = y2 + g2 .* ny2;

        x3geom = x3 + g3 .* nx3;
        y3geom = y3 + g3 .* ny3;

        hold on
        for i = 1 : length(x1)
            plot([x1(i), x1geom(i)], [y1(i), y1geom(i)], 'c-')
            plot([x2(i), x2geom(i)], [y2(i), y2geom(i)], 'c-')
            plot([x3(i), x3geom(i)], [y3(i), y3geom(i)], 'c-')
        end % for i
        hold off
    end % if

    % compute trimming change at node A (using end of edge 3 and beg of edge 1)
    dt   = [        ny3(end),     -ny1(1);       nx3(end),       -nx1(1  )] \ ...
           [g3(end)*nx3(end)-g1(1)*nx1(1); g1(1)*ny1(1  )-g3(end)*ny3(end)];
    dt3A = dt(1);
    dt1A = dt(2);

    % compute trimming change at node B (using end of edge 1 and beg of edge 2)
    dt   = [        ny1(end),     -ny2(1);       nx1(end),       -nx2(1  )] \ ...
           [g1(end)*nx1(end)-g2(1)*nx2(1); g2(1)*ny2(1  )-g1(end)*ny1(end)];
    dt1B = dt(1);
    dt2B = dt(2);

    % compute trimming change at node C (using end of edge 2 and beg of edge 3)
    dt   = [        ny2(end),     -ny3(1);       nx2(end),       -nx3(1  )] \ ...
           [g2(end)*nx2(end)-g3(1)*nx3(1); g3(1)*ny3(1  )-g2(end)*ny2(end)];
    dt2C = dt(1);
    dt3C = dt(2);

    % compute and show the tessellation sensitivities
    if (iplot >= 2)
        tt     = linspace(0, 1, length(x1));

        x1tess = x1geom - ((1-tt) * dt1A + tt * dt1B) .* ny1;
        y1tess = y1geom + ((1-tt) * dt1A + tt * dt1B) .* nx1;

        x2tess = x2geom - ((1-tt) * dt2B + tt * dt2C) .* ny2;
        y2tess = y2geom + ((1-tt) * dt2B + tt * dt2C) .* nx2;

        x3tess = x3geom - ((1-tt) * dt3C + tt * dt3A) .* ny3;
        y3tess = y3geom + ((1-tt) * dt3C + tt * dt3A) .* nx3;

        hold on
        for i = 1 : length(x1)
            plot([x1(i), x1tess(i)], [y1(i), y1tess(i)], 'm-')
            plot([x2(i), x2tess(i)], [y2(i), y2tess(i)], 'm-')
            plot([x3(i), x3tess(i)], [y3(i), y3tess(i)], 'm-')
        end % for i
        hold off
    end % if

end % function designV1
