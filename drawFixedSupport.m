% ------------------------------------------------------------
% drawFixedSupport
% ------------------------------------------------------------
% PURPOSE:
%   Draws a small graphical symbol representing a fixed (built‑in)
%   support at coordinates (x, y) on a 2‑D frame or beam plot.
%
% USAGE:
%   drawFixedSupport(x, y, span)
%
% INPUTS:
%   x, y   : Coordinates of the base point of the support
%   span   : Reference length (e.g., typical member span) used to
%            scale the size of the support symbol
%
% DESCRIPTION:
%   ‑  The function draws a vertical rectangular block beneath the
%      node location and fills it with diagonal hatch lines to
%      represent a fully fixed support.
%   ‑  The patch dimensions are proportional to 'span' so the symbol
%      scales automatically with the frame size.
%
% NOTES:
%   ‑  Purely graphical; does not influence analysis results.
%   ‑  Designed to be called after plotting members and joints, e.g.:
% ------------------------------------------------------------

function drawFixedSupport(x,y,span)
    w = 0.04*span; 
    h = 0.02*span;

    % Single solid block under the joint
    patch([x-w x+w x+w x-w], [y y y-h y-h], ...
        [0.65 0.65 0.65], 'EdgeColor','k', 'LineWidth',1.2);
end