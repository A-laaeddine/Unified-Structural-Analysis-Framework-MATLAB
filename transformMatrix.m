% ------------------------------------------------------------
% transformMatrix
% ------------------------------------------------------------
% PURPOSE:
%   Builds the 6×6 transformation matrix for a 2‑D frame element
%   used to convert vectors or matrices between the local and
%   global coordinate systems.
%
% USAGE:
%   T = transformMatrix(c, s)
%
% INPUTS:
%   c  : Cosine of the element inclination angle (cosθ)
%   s  : Sine of the element inclination angle (sinθ)
%
% OUTPUT:
%   T  : 6×6 transformation matrix such that
%        [k_global] = [T]' * [k_local] * [T]
%        and {F_global} = [T]' * {F_local}, etc.
%
% DESCRIPTION:
%   - Constructs the 3×3 rotation sub‑matrix for a 2‑D element:
%         t = [ c  s  0;
%              -s  c  0;
%               0  0  1 ];
%   - Uses blkdiag(t, t) to place this rotation in both end‑node
%     blocks, producing a full 6×6 element transformation matrix.
%
% NOTES:
%   - Sign convention assumes global Y is upward and element
%     local x‑axis runs from start node to end node.
%   - For in‑plane frame analysis only (2 translations + 1 rotation
%     per node).
% ------------------------------------------------------------

function T = transformMatrix(c, s)
    t = [c s 0; -s c 0; 0 0 1];
    T = blkdiag(t, t);
end