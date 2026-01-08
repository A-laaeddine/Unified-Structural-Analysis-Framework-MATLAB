    % ------------------------------------------------------------
    % localStiffness
    % ------------------------------------------------------------
    % PURPOSE:
    %   Computes the 6×6 local stiffness matrix of a 2-D frame (beam-column)
    %   element accounting for axial and flexural behavior.
    %
    % USAGE:
    %   k = localStiffness(E, I, A, L)
    %
    % INPUTS:
    %   E : Young's modulus of the member material
    %   I : Second moment of area about the local bending axis
    %   A : Cross-sectional area
    %   L : Member length
    %
    % OUTPUT:
    %   k : 6×6 local element stiffness matrix in the member's
    %       local coordinate system, ordered as:
    %       [u1 v1 θ1 u2 v2 θ2]
    %
    % DESCRIPTION:
    %   - Includes axial stiffness (AE/L)
    %   - Assumes small deformations and linear elastic behavior.
    %   - The matrix relates local nodal forces to local nodal
    %     displacements before any coordinate transformation.
    %
    % NOTES:
    %   - Valid for prismatic members with constant E, A, and I.
    %   - Shear deformation effects are neglected (classical beam theory).
    %   - Typically transformed to the global coordinate system
    %     prior to assembly into the structure stiffness matrix.
    % ------------------------------------------------------------
    
    function k = localStiffness(E, I, A, L)
        k = [ A*E/L 0 0 -A*E/L 0 0; 
              0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
              0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
              -A*E/L 0 0 A*E/L 0 0;
              0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
              0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L ];
    end