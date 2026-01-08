% ------------------------------------------------------------
% manualInput
% ------------------------------------------------------------
% PURPOSE:
%   Provides an interactive console-based method for defining all
%   geometric, material, loading, and boundary-condition data for
%   a 2-D frame analysis model.
%
% USAGE:
%   [joints, members, jLoads, mLoads, supports] = manualInput()
%
% OUTPUTS:
%   joints   : [ID, x, y] joint coordinates
%   members  : [ID, startJoint, endJoint, E, I, A] member connectivity
%              and section/material properties
%   jLoads   : [jointID, Fx, Fy, M] applied joint forces and moments
%   mLoads   : [memberID, wy] uniformly distributed member loads
%   supports : [jointID, Ux, Uy, Rot] support constraints
%              (1 = restrained, 0 = free)
%
% DESCRIPTION:
%   - Prompts the user to manually enter each data set via the MATLAB
%     command window.
%   - Automatically assigns predefined default examples when the
%     user provides empty input, enabling rapid testing and debugging.
%   - Ensures consistent data structures for downstream assembly,
%     stiffness formulation, and visualization routines.
%
% NOTES:
%   - Intended primarily for small problems, demonstrations,
%     or educational use.
%   - Performs no validation beyond empty-input handling.
% ------------------------------------------------------------

function [joints, members, jLoads, mLoads, supports] = manualInput()

% 1. JOINTS
% Format: [ID, x, y]
joints = input('Joints [ID x y; ...]: ');
if isempty(joints)
    error('Input Error: The "joints" matrix cannot be empty. The structure must have at least two joints.');
end

% 2. MEMBERS
% Format: [ID, Start, End, E, I, A]
members = input('Members [ID Start End E I A; ...]: ');
if isempty(members)
    error('Input Error: The "members" matrix cannot be empty. A structure requires at least one member.');
end

% 3. JOINT LOADS
% Format: [JointID, Fx, Fy, M]
% (Keeping empty allowed here, as a structure can exist without joint loads)
jLoads = input('Joint Loads [JointID Fx Fy M; ...]: ');

% 4. MEMBER LOADS
% Format: [MemberID, LoadValue, Type, Dir]
% (Keeping empty allowed here)
mLoads = input('Member Loads [MemberID val type dir; ...]: '); % for dir: 1-> Perpendicular to local member axis
                                                               %          2-> Along member local axis
                                                                        
% 5. SUPPORTS
% Format: [JointID, Ux, Uy, Rot] (1=fixed, 0=free)
supports = input('Supports [JointID Ux Uy Rot; ...]: ');
if isempty(supports)
    error('Stability Error: The "supports" matrix is empty. The structure is unstable (unconstrained) and cannot be solved.');
end

end
