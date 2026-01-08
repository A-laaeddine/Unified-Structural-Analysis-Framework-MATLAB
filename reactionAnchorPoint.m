% ------------------------------------------------------------
% reactionAnchorPoint
% ------------------------------------------------------------
% PURPOSE:
%   Determines an appropriate starting point for plotting support
%   reaction arrows and labels so they do not overlap with the
%   graphical support symbol.
%
% USAGE:
%   [xr, yr] = reactionAnchorPoint(x, y, span, fixX, fixY, fixR)
%
% INPUTS:
%   x, y        : Coordinates of the joint location
%   span       : Reference length used to scale positional offsets
%   fixX       : X-direction restraint flag (1 = restrained, 0 = free)
%   fixY       : Y-direction restraint flag (1 = restrained, 0 = free)
%   fixR       : Rotational restraint flag (1 = restrained, 0 = free)
%
% OUTPUTS:
%   xr, yr     : Coordinates of the recommended anchor point for
%                reaction arrows and text labels
%
% DESCRIPTION:
%   - Selects a vertical offset based on the support type
%     (fixed, pin, roller, or generic).
%   - Places reaction graphics slightly below the support symbol
%     to maintain visual clarity and avoid overlap.
%   - Uses span-proportional offsets to ensure consistent spacing
%     across different model scales.
%
% NOTES:
%   - Purely graphical; does not affect numerical results.
%   - Intended for use in post-processing and visualization
%     routines that display support reactions.
% ------------------------------------------------------------

function [xr, yr] = reactionAnchorPoint(x, y, span, fixX, fixY, fixR)
    if fixX==1 && fixY==1 && fixR==1
        % fixed support block extends below joint => place reactions below it
        yr = y - 0.13*span;
        xr = x;
    elseif fixX==1 && fixY==1 && fixR==0
        % pin => place slightly below
        yr = y - 0.08*span;
        xr = x;
    elseif (fixX==0 && fixY==1 && fixR==0) || (fixX==1 && fixY==0 && fixR==0)
        % roller => place slightly below
        yr = y - 0.08*span;
        xr = x;
    else
        % generic
        yr = y - 0.06*span;
        xr = x;
    end
end