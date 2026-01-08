% ------------------------------------------------------------
% drawPinSupport
% ------------------------------------------------------------
% PURPOSE:
%   Draws a graphical symbol representing a pin (hinged) support
%   at a specified joint location in a 2-D frame or beam plot.
%
% USAGE:
%   drawPinSupport(x, y, span)
%
% INPUTS:
%   x, y : Coordinates of the support location
%   span : Reference length used to scale the symbol size
%
% DESCRIPTION:
%   - Renders a triangular support block beneath the joint to
%     represent translational restraint with rotational freedom.
%   - Draws a short base line under the triangle to indicate
%     contact with the supporting surface.
%   - All geometric proportions are scaled relative to the overall
%     frame span for consistent appearance.
%
% NOTES:
%   - Graphical utility only; does not influence analysis results.
%   - Designed to complement fixed and roller support symbols
%     within frame visualization routines.
% ------------------------------------------------------------

function drawPinSupport(x,y,span)
    s = 0.05*span;
    patch([x x-s x+s],[y-0.02*span y-0.10*span y-0.10*span], ...
        [0.85 0.85 0.85],'EdgeColor','k','LineWidth',1.2);
    plot([x-1.2*s x+1.2*s],[y-0.10*span y-0.10*span],'k-','LineWidth',1);
end