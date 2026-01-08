% =========================================================================
% 2D STRUCTURAL ANALYSIS SOLVER
% (DIRECT STIFFNESS METHOD – BEAMS, FRAMES & TRUSSES)
% Developed by: Mohamad Alaaeddine
% =========================================================================
% DESCRIPTION:
% This program performs linear static analysis of 2D structural systems,
% including beams, frames, and trusses, using the Direct Stiffness Method.
% Element stiffness matrices are formulated in local coordinates and
% transformed into the global coordinate system prior to assembly.
%
% The solver supports two analysis modes:
% 1. Beam / Frame analysis (axial + bending + rotation DOFs).
% 2. Truss analysis (axial DOFs only).
%
% The program computes:
% - Joint displacements,
% - Support reactions,
% - (For trusses) internal axial member forces.
%
% The solver accounts for:
% 1. Concentrated joint loads (forces and moments).
% 2. Distributed member loads converted to equivalent joint loads (Pf).
% 3. Support restraints (fixed, pinned, or roller boundary conditions).
%
% -------------------------------------------------------------------------
% IMPORTANT NOTES & CONVENTIONS
% -------------------------------------------------------------------------
% 1. Unit Consistency:
%    All quantities must be defined in a single coherent unit system.
%    Typical choices include:
%    • Kip–ft system
%    • kN–m system
%
% 2. Global Coordinate System:
%    Global X → horizontal to the right
%    Global Y → vertical upward
%
% 3. Load & Support Sign Convention:
%    • Positive Fx acts to the right.
%    • Positive Fy acts upward.
%    • Positive Mz (moment) is counter-clockwise.
%    • In the supports matrix: 1 = restrained, 0 = free.
%    • For pinned or roller supports, set Rot = 0.
%
% -------------------------------------------------------------------------
% LIMITATIONS & MODELING GUIDELINES
% -------------------------------------------------------------------------
% • Applicable only to 2D structural systems.
% • Linear elastic, small-displacement, static analysis only.
% • Beam and frame member loading is limited to:
%     – uniformly distributed loads (UDL),
%     – concentrated member loads applied at mid-span.
% • Concentrated member loads not acting at mid-span must be modeled
%   by subdividing the member and applying the load as a joint force.
% • Inclined UDL along members must be divided into x- and y- components
%   before inputting them into the program.
% • Assumes constant EI along each member. For pure beam problems,
%   E = 1 and I = 1 may be used as normalized values.
% • For beam-dominated behavior, a large cross-sectional area A should
%   be used to minimize axial deformation effects.
% • Truss elements carry axial forces only; bending effects are ignored.
% • Nonlinear material behavior, geometric nonlinearity, stability
%   (buckling), dynamic effects, and time-dependent loads are not included.
% =========================================================================

%% 0. CLEAR & HEADER
clc; clear; close all;

fprintf('===============================================\n');
fprintf('     FRAME ANALYSIS SOLVER (DIRECT STIFFNESS)  \n');
fprintf('           Made by Mohamad Alaaeddine          \n');
fprintf('===============================================\n');

%% 1. STRUCTURE TYPE SELECTION
fprintf('Select Structure Type:\n');
fprintf('1. Frame or Beam\n');
fprintf('2. Truss\n');
structureType = input('Enter choice (1 or 2): ');

%% ===============================================================
%                    FRAME / BEAM ANALYSIS
% ===============================================================
if structureType == 1

fprintf('\nFrame or Beam analysis mode selected.\n');

%% 2. DATA INPUT
[joints, members, jLoads, mLoads, supports] = manualInput();

%% 3. INITIALIZATION
nJ   = size(joints, 1);
nM   = size(members, 1);
ndof = 3 * nJ;

S       = zeros(ndof);      % Global stiffness matrix
F_joint = zeros(ndof, 1);   % Joint load vector (P)
Pf      = zeros(ndof, 1);   % Equivalent member load vector (Pf)

%% 4. GLOBAL STIFFNESS ASSEMBLY & MEMBER LOADS
for e = 1:nM

    % --- Member Identification ---
    id = members(e,1);
    n1 = members(e,2);
    n2 = members(e,3);
    E  = members(e,4);
    I  = members(e,5);
    A  = members(e,6);

    % --- Geometry ---
    x1 = joints(n1,2);  y1 = joints(n1,3);
    x2 = joints(n2,2);  y2 = joints(n2,3);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2-x1)/L;
    s = (y2-y1)/L;

    % --- Local & Transformation Matrices ---
    kL = localStiffness(E, I, A, L);
    T  = transformMatrix(c, s);

    dofs = [3*n1-2, 3*n1-1, 3*n1, ...
            3*n2-2, 3*n2-1, 3*n2];

    % --- Assemble Global Stiffness ---
    S(dofs,dofs) = S(dofs,dofs) + T' * kL * T;

    % --- Member Load Logic ---
    Qf = zeros(6,1);
    idx = find(mLoads(:,1) == id, 1);

    if ~isempty(idx) && size(mLoads,2) >= 4
        val  = mLoads(idx,2);
        type = mLoads(idx,3);
        dir  = mLoads(idx,4);

        if dir == 1
            p_trans = val * s;
            p_axial = val * c;
        else
            p_trans = -val * c;
            p_axial =  val * s;
        end

        if type == 1
            Qf = [ p_axial*L/2;
                   p_trans*L/2;
                   p_trans*L^2/12;
                   p_axial*L/2;
                   p_trans*L/2;
                  -p_trans*L^2/12 ];
        elseif type == 2
            Qf = [ p_axial/2;
                   p_trans/2;
                   p_trans*L/8;
                   p_axial/2;
                   p_trans/2;
                  -p_trans*L/8 ];
        end
    end

    Pf(dofs) = Pf(dofs) + T' * Qf;
end

%% 5. JOINT LOAD VECTOR
for i = 1:size(jLoads,1)
    jID = jLoads(i,1);
    F_joint(3*jID-2:3*jID) = ...
        F_joint(3*jID-2:3*jID) + jLoads(i,2:4)';
end

%% 6. BOUNDARY CONDITIONS & SOLUTION
fixed = [];
for i = 1:size(supports,1)
    sn = supports(i,1);
    if supports(i,2)==1, fixed = [fixed, 3*sn-2]; end
    if supports(i,3)==1, fixed = [fixed, 3*sn-1]; end
    if supports(i,4)==1, fixed = [fixed, 3*sn];   end
end

free = setdiff(1:ndof, fixed);

d = zeros(ndof,1);
if ~isempty(free)
    d(free) = S(free,free) \ (F_joint(free) - Pf(free));
end

%% 7. REACTIONS
R = S*d + Pf - F_joint;

%% 8. OUTPUT RESULTS
fprintf('\n===============================================\n');
fprintf('                ANALYSIS RESULTS               \n');
fprintf('===============================================\n');

fprintf('\nJOINT DISPLACEMENTS:\n');
fprintf('Joint |      Ux      |      Uy      |     Rot      \n');
for i = 1:nJ
    fprintf('%5d | %12.4e | %12.4e | %12.4e\n', ...
        i, d(3*i-2), d(3*i-1), d(3*i));
end

fprintf('\nSUPPORT REACTIONS:\n');
fprintf('Joint |      Fx      |      Fy      |      Mz      \n');
for i = 1:size(supports,1)
    j = supports(i,1);
    fprintf('%5d | %12.4f | %12.4f | %12.4f\n', ...
        j, R(3*j-2), R(3*j-1), R(3*j));
end

%% 9. PLOTTING
plotFrameAndForces(joints, members, jLoads, mLoads, supports, R);

%% ===============================================================
%                         TRUSS ANALYSIS
% ===============================================================
elseif structureType == 2

fprintf('\nTruss analysis mode selected.\n');

%% 2. DATA INPUT
[joints, members, jLoads, mLoads, supports] = manualInput();

%% 3. INITIALIZATION
nJ   = size(joints,1);
nM   = size(members,1);
ndof = 2 * nJ;

S = zeros(ndof);
P = zeros(ndof,1);

%% 4. GLOBAL STIFFNESS ASSEMBLY
for e = 1:nM

    id = members(e,1);
    n1 = members(e,2);
    n2 = members(e,3);
    E  = members(e,4);
    A  = members(e,6);

    x1 = joints(n1,2); y1 = joints(n1,3);
    x2 = joints(n2,2); y2 = joints(n2,3);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2-x1)/L;
    s = (y2-y1)/L;

    kL = (E*A/L) * [1 -1; -1 1];
    T  = [c s 0 0; 0 0 c s];

    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    S(dofs,dofs) = S(dofs,dofs) + T' * kL * T;
end

%% 5. LOAD VECTOR
for i = 1:size(jLoads,1)
    jID = jLoads(i,1);
    P(2*jID-1:2*jID) = P(2*jID-1:2*jID) + jLoads(i,2:3)';
end

%% 6. BOUNDARY CONDITIONS & SOLUTION
fixed = [];
for i = 1:size(supports,1)
    sn = supports(i,1);
    if supports(i,2)==1, fixed = [fixed, 2*sn-1]; end
    if supports(i,3)==1, fixed = [fixed, 2*sn];   end
end

free = setdiff(1:ndof, fixed);

d = zeros(ndof,1);
if ~isempty(free)
    d(free) = S(free,free) \ P(free);
end

%% 7. MEMBER FORCES & REACTIONS
R = S*d - P;

fprintf('\n===============================================\n');
fprintf('                ANALYSIS RESULTS               \n');
fprintf('===============================================\n');

fprintf('\nINTERNAL MEMBER FORCES (Axial):\n');
fprintf('Member |    Force (k or kN)    |   Type   \n');

for e = 1:nM
    n1 = members(e,2);
    n2 = members(e,3);
    E  = members(e,4);
    A  = members(e,6);

    L = sqrt((joints(n2,2)-joints(n1,2))^2 + ...
             (joints(n2,3)-joints(n1,3))^2);

    c = (joints(n2,2)-joints(n1,2))/L;
    s = (joints(n2,3)-joints(n1,3))/L;

    u_global = [d(2*n1-1); d(2*n1); d(2*n2-1); d(2*n2)];

    axialForce = (E*A/L) * [-c -s c s] * u_global;

    type = 'Tension';
    if axialForce < 0, type = 'Compression'; end

    fprintf('%6d | %20.4f | %s\n', e, abs(axialForce), type);
end

fprintf('\nSUPPORT REACTIONS (Truss):\n');
fprintf('Joint |      Fx      |      Fy      \n');

for i = 1:size(supports,1)
    jID = supports(i,1);
    fprintf('%5d | %12.4f | %12.4f\n', ...
        jID, R(2*jID-1), R(2*jID));
end

%% 8. PLOTTING
plotTrussAndForces(joints, members, jLoads, supports, R);

end