% ------------------------------------------------------------
% drawMomentArc
% ------------------------------------------------------------
% PURPOSE:
%   Draws a curved arrow (moment symbol) at a specified location
%   on a 2‑D frame plot.  The arrow direction indicates the
%   sign of the bending moment: counter‑clockwise (CCW) for
%   positive values, clockwise (CW) for negative values.
%
% USAGE:
%   drawMomentArc(x, y, r, M, color)
%
% INPUTS:
%   x, y   : Coordinates of the arc center (location of moment)
%   r      : Radius of the moment‑arc symbol
%   M      : Moment magnitude (positive → CCW, negative → CW)
%   color  : RGB triplet or MATLAB color character (e.g. 'r','k')
%
% DESCRIPTION:
%   ‑  The function generates an arc segment (using parametric
%      angles theta) centered at (x,y).
%   ‑  The arc extent and orientation follow the sign of M.
%   ‑  A small filled triangular arrowhead is attached at the
%      end of the arc to indicate moment direction.
%
% NOTES:
%   ‑  Intended for plotting schematic moment symbols on frame
%      geometry or deformed‑shape figures.
%   ‑  This routine is purely graphical and does not scale
%      arrow size by actual moment value.
% ------------------------------------------------------------

function drawMomentArc(x,y,r,M,color)
    if M > 0
        th = linspace(0.2*pi, 1.6*pi, 80);
    else
        th = linspace(1.8*pi, 0.4*pi, 80);
    end
    xp = x + r*cos(th); yp = y + r*sin(th);
    plot(xp, yp, 'Color', color, 'LineWidth',2);

    xe=xp(end); ye=yp(end);
    ang = atan2(yp(end)-yp(end-1), xp(end)-xp(end-1));
    ah = 0.25*r;
    p1=[xe ye];
    p2=[xe-ah*cos(ang-pi/6) ye-ah*sin(ang-pi/6)];
    p3=[xe-ah*cos(ang+pi/6) ye-ah*sin(ang+pi/6)];
    fill([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)], color, 'EdgeColor', color);
end