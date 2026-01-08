% ------------------------------------------------------------
% plotEndForceSmall
% ------------------------------------------------------------
% PURPOSE:
%   Draws a small force vector at a specified point to represent
%   an element end force in a 2-D frame or beam visualization.
%
% USAGE:
%   plotEndForceSmall(x, y, q, span)
%
% INPUTS:
%   x, y : Coordinates of the force application point
%   q    : Force vector [Fx, Fy] in global coordinates
%   span : Reference length used to scale the arrow size
%
% DESCRIPTION:
%   - Computes the force magnitude and direction from the input
%     vector components.
%   - Plots a fixed-length arrow oriented along the force direction,
%     ensuring visual consistency across the structure.
%   - Intended for compact representation of internal or end forces
%     without dominating the overall plot.
%
% NOTES:
%   - Graphical utility only; does not affect analysis results.
%   - Forces with zero magnitude are ignored.
%   - Typically used alongside member or joint force diagrams
%     for clarity in post-processing.
% ------------------------------------------------------------

function plotEndForceSmall(x,y,q,span)
    Fx=q(1); Fy=q(2);
    Fmag=hypot(Fx,Fy);
    if Fmag>0
        ux=Fx/Fmag; uy=Fy/Fmag;
        L=0.10*span;
        quiver(x,y,ux*L,uy*L,0,'Color',[0.1 0.6 0.2],'LineWidth',2,'MaxHeadSize',0.6);
    end
end