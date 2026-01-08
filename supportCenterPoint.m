% ------------------------------------------------------------
% supportCenterPoint
% ------------------------------------------------------------
% PURPOSE:
%   Computes the geometric center of a support symbol so that
%   reaction forces and moments are drawn in a visually centered
%   and physically intuitive location.
%
% USAGE:
%   [xc, yc] = supportCenterPoint(x, y, span, fixX, fixY, fixR, fixedH)
%
% INPUTS:
%   x, y     : Coordinates of the joint location
%   span     : Reference length used to scale positional offsets
%   fixX     : X-direction restraint flag (1 = restrained, 0 = free)
%   fixY     : Y-direction restraint flag (1 = restrained, 0 = free)
%   fixR     : Rotational restraint flag (1 = restrained, 0 = free)
%   fixedH   : Height of the fixed support block (used only for
%              fully fixed supports)
%
% OUTPUTS:
%   xc, yc   : Coordinates of the geometric center of the support symbol
%
% DESCRIPTION:
%   - Identifies the support type based on restraint flags.
%   - For fixed supports, returns the midpoint of the rectangular
%     support block.
%   - For pin and roller supports, returns a point slightly below
%     the joint to align reactions with the symbol center.
%   - Uses span-proportional offsets to ensure consistent placement
%     across different model scales.
%
% NOTES:
%   - Graphical utility only; does not influence analysis results.
%   - Intended for post-processing routines that plot reactions
%     and support moments without visual overlap.
% ------------------------------------------------------------

function [xc, yc] = supportCenterPoint(x, y, span, fixX, fixY, fixR, fixedH)
    if fixX==1 && fixY==1 && fixR==1
        % fixed block is drawn from y down to y-fixedH
        xc = x;
        yc = y - fixedH/2;
    elseif fixX==1 && fixY==1 && fixR==0
        % pin triangle is below joint: center slightly below
        xc = x;
        yc = y - 0.075*span;
    elseif (fixX==0 && fixY==1 && fixR==0) || (fixX==1 && fixY==0 && fixR==0)
        % roller symbol is below joint: center slightly below
        xc = x;
        yc = y - 0.075*span;
    else
        % generic
        xc = x;
        yc = y - 0.06*span;
    end
end