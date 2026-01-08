% ------------------------------------------------------------
% drawRollerSupport
% ------------------------------------------------------------
% PURPOSE:
%   Draws a graphical symbol representing a roller support at a
%   specified joint location in a 2-D frame or beam plot.
%
% USAGE:
%   drawRollerSupport(x, y, span)
%
% INPUTS:
%   x, y : Coordinates of the support location
%   span : Reference length used to scale the symbol size
%
% DESCRIPTION:
%   - Renders a short horizontal base line beneath the joint.
%   - Adds two circular rollers below the base to indicate
%     translational restraint in one direction and freedom in the
%     orthogonal direction.
%   - All dimensions are proportional to the global frame span to
%     maintain visual consistency across different model scales.
%
% NOTES:
%   - Purely graphical; does not affect structural analysis results.
%   - Intended to be used alongside pin and fixed support symbols
%     for clear boundary-condition visualization.
% ------------------------------------------------------------

function drawRollerSupport(x,y,span)
    s = 0.05*span;
    plot([x-s x+s],[y-0.08*span y-0.08*span],'k-','LineWidth',1.2);
    plot(x-0.5*s, y-0.10*span,'ko','MarkerFaceColor','w','MarkerSize',6);
    plot(x+0.5*s, y-0.10*span,'ko','MarkerFaceColor','w','MarkerSize',6);
end