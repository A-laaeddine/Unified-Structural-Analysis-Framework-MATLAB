% ------------------------------------------------------------
% plotTrussAndForces
% ------------------------------------------------------------
% PURPOSE:
%   Plots the complete 2-D truss geometry together with supports,
%   applied joint forces (Fx, Fy), and computed support reactions
%   (Rx, Ry) in a clean, scaled visualization.
%
% USAGE:
%   plotTrussAndForces(joints, members, jLoads, supports, R)
%
% INPUTS:
%   joints   : [ID, x, y] table of joint identifiers and coordinates
%   members  : [ID, joint_i, joint_j, ...] truss member connectivity
%   jLoads   : [jointID, Fx, Fy] applied joint forces
%   supports : [jointID, fixX, fixY] support conditions
%              (1 = restrained, 0 = free)
%   R        : global reaction vector (2*nJoints × 1) from truss solver,
%              ordered as [Rx1, Ry1, Rx2, Ry2, ...]
%
% DESCRIPTION:
%   - Draws all truss members as straight axial elements between joints.
%   - Displays joints as nodes with automatic scaling based on structure span.
%   - Renders pinned and roller supports according to restraint conditions.
%   - Plots joint force components (Fx and Fy) as arrows originating at joints.
%   - Plots support reactions (Rx and Ry) at restrained joints only.
%   - Adds stacked numerical labels for applied loads and reactions for
%     clear interpretation of magnitude and direction.
%
% NOTES:
%   - Visualization utility only; does not affect analysis results.
%   - Intended strictly for 2-D truss systems (2 DOF per joint: Fx, Fy).
%   - Robust to missing or invalid joint/member references
%     (warnings are issued and plotting continues).
%   - Arrow lengths are automatically scaled for consistent appearance
%     regardless of load magnitude.
% ------------------------------------------------------------

function plotTrussAndForces(joints, members, jLoads, supports, R)

tol = 1e-12;

% 0. INPUT CHECKS
if isempty(joints) || isempty(members)
    error('JOINTS and MEMBERS must not be empty.');
end

if numel(unique(joints(:,1))) ~= size(joints,1)
    error('JOINTS has duplicate Joint IDs.');
end

% 1. JOINT ID → ROW INDEX MAP
id2row = containers.Map(num2cell(joints(:,1)), num2cell(1:size(joints,1)));

% 2. FIGURE SETUP
figure('Color','w', 'Name','Truss Analysis', 'NumberTitle','off');
clf; hold on; axis equal;
axis off; grid off;
title('Truss Analysis');

xs = joints(:,2); ys = joints(:,3);
span = max(max(xs)-min(xs), max(ys)-min(ys));
if span == 0, span = 1; end
margin = 0.30*span;
xlim([min(xs)-margin, max(xs)+margin]);
ylim([min(ys)-margin, max(ys)+margin]);

% ===============================================================
% 3. PLOT MEMBERS & JOINTS
% ===============================================================
for e = 1:size(members,1)
    j1id = members(e,2);
    j2id = members(e,3);

    if ~isKey(id2row,j1id) || ~isKey(id2row,j2id)
        warning('Member %g references missing joint ID(s): %g, %g', members(e,1), j1id, j2id);
        continue;
    end

    r1 = id2row(j1id);
    r2 = id2row(j2id);

    x1 = joints(r1,2); y1 = joints(r1,3);
    x2 = joints(r2,2); y2 = joints(r2,3);

    plot([x1 x2],[y1 y2], '-', 'Color',[0.1 0.4 0.9], 'LineWidth',2);
end

plot(xs, ys, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);

% ===============================================================
% 4. PLOT SUPPORTS (PIN / ROLLER)
% ===============================================================
if ~isempty(supports)
    for i=1:size(supports,1)
        jid = supports(i,1);
        if ~isKey(id2row,jid), continue; end
        rr = id2row(jid);

        x = joints(rr,2);
        y = joints(rr,3);

        fixX = supports(i,2);
        fixY = supports(i,3);

        if fixX==1 && fixY==1
            drawPinSupport(x,y,span);
        elseif (fixX==0 && fixY==1) || (fixX==1 && fixY==0)
            drawRollerSupport(x,y,span);
        else
            plot(x,y,'ks','MarkerFaceColor',[0.75 0.75 0.75],'MarkerSize',7);
        end
    end
end

% ===============================================================
% 5. PLOT JOINT LOADS
% ===============================================================
if ~isempty(jLoads)
    % Scaling references for Fx and Fy separately
    FxAbs = abs(jLoads(:,2)); FxPos = FxAbs(FxAbs>tol);
    FyAbs = abs(jLoads(:,3)); FyPos = FyAbs(FyAbs>tol);

    FxRef = 1; FyRef = 1;
    if ~isempty(FxPos)
        FxRef = prctile(FxPos,75);
        if FxRef==0, FxRef=max(FxPos); end
        if FxRef==0, FxRef=1; end
    end
    if ~isempty(FyPos)
        FyRef = prctile(FyPos,75);
        if FyRef==0, FyRef=max(FyPos); end
        if FyRef==0, FyRef=1; end
    end

    Lmax = 0.18*span;
    Lmin = 0.05*span;

    for i=1:size(jLoads,1)
        jid = jLoads(i,1);
        if ~isKey(id2row,jid)
            warning('J-Load references missing joint ID: %g', jid);
            continue;
        end
        rr = id2row(jid);

        x = joints(rr,2);
        y = joints(rr,3);

        Fx = jLoads(i,2);
        Fy = jLoads(i,3);

        % Fx arrow (horizontal)
        if abs(Fx) > tol
            Lx = Lmax*(abs(Fx)/(abs(Fx)+FxRef));
            Lx = max(Lx,Lmin);
            sgn = sign(Fx);
            quiver(x, y, sgn*Lx, 0, 0, 'Color',[0.95 0.75 0.1], ...
                   'LineWidth',2, 'MaxHeadSize',0.6);
        end

        % Fy arrow (vertical)
        if abs(Fy) > tol
            Ly = Lmax*(abs(Fy)/(abs(Fy)+FyRef));
            Ly = max(Ly,Lmin);
            sgn = sign(Fy);
            quiver(x, y, 0, sgn*Ly, 0, 'Color',[0.95 0.75 0.1], ...
                   'LineWidth',2, 'MaxHeadSize',0.6);
        end

        % stacked label under joint
        if (abs(Fx) > tol) || (abs(Fy) > tol)
            txtL = sprintf('Fx = %.3g\nFy = %.3g', Fx, Fy);
            text(x, y - 0.12*span, txtL, ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', ...
                'FontSize',9, 'FontWeight','bold', ...
                'Color',[0.55 0.42 0.02], ...
                'BackgroundColor','none', 'EdgeColor','none', ...
                'Margin',1, 'Clipping','on');
        end
    end
end

% ===============================================================
% 6. PLOT REACTIONS
% ===============================================================
if ~isempty(R) && ~isempty(supports)
    rxLen = 0.10*span;
    ryLen = 0.10*span;

    for i = 1:size(supports,1)
        jid = supports(i,1);
        if ~isKey(id2row,jid), continue; end
        rr = id2row(jid);

        x = joints(rr,2);
        y = joints(rr,3);

        fixX = supports(i,2);
        fixY = supports(i,3);

        % Rx and Ry components (2 DOF per joint)
        Rx = R(2*rr-1);
        Ry = R(2*rr);

        % Rx arrow
        if fixX == 1 && abs(Rx) > tol
            sgn = sign(Rx);
            quiver(x, y, sgn*rxLen, 0, 0, 'Color',[0.85 0.1 0.1], ...
                   'LineWidth',2.2, 'MaxHeadSize',0.7);
        end

        % Ry arrow
        if fixY == 1 && abs(Ry) > tol
            sgn = sign(Ry);
            quiver(x, y, 0, sgn*ryLen, 0, 'Color',[0.85 0.1 0.1], ...
                   'LineWidth',2.2, 'MaxHeadSize',0.7);
        end

        % stacked reaction label under support
        RxDisp = Rx; if fixX~=1, RxDisp=0; end
        RyDisp = Ry; if fixY~=1, RyDisp=0; end

        txtR = sprintf('Rx = %.3g\nRy = %.3g', RxDisp, RyDisp);
        text(x, y - 0.16*span, txtR, ...
             'HorizontalAlignment','center', 'VerticalAlignment','top', ...
             'FontSize',9, 'FontWeight','bold', ...
             'Color',[0.55 0 0], ...
             'BackgroundColor','none', 'EdgeColor','none', ...
             'Margin',1, 'Clipping','on');
    end
end

hold off;
end
