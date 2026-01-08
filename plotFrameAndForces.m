% ------------------------------------------------------------
% plotFrameAndForces
% ------------------------------------------------------------
% PURPOSE:
%   Plots the complete 2-D frame geometry together with supports,
%   joint loads, joint moments, member distributed load markers,
%   AND support reactions (Rx, Ry, M) in a clean, scaled visualization.
%
% USAGE:
%   plotFrameAndForces(joints, members, jLoads, mLoads, supports, R)
%
% INPUTS:
%   joints   : [ID, x, y] table of joint identifiers and coordinates
%   members  : [ID, joint_i, joint_j, ...] connectivity of frame members
%   jLoads   : [jointID, Fx, Fy, Mz] joint force and moment loads
%   mLoads   : [memberID, val, type, ldir]  type: 1=UDL, 2=Point ; ldir: 1=X, 2=Y
%   supports : [jointID, fixX, fixY, fixR] support conditions (1 restrained, 0 free)
%   R        : global reaction vector (ndof x 1) from solver: R = S*d + Pf - F_nodal
%
% NOTES:
%   - Visualization only; robust to missing IDs (warnings issued, plotting continues).
%   - IMPORTANT: Reaction indexing assumes Joint IDs are 1..n (same as your solver DOF numbering).
% ------------------------------------------------------------

function plotFrameAndForces(joints, members, jLoads, mLoads, supports, R)
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
figure('Color','w', 'Name','Structure Analysis', 'NumberTitle','off');
clf; hold on; axis equal;
axis off; grid off;
title('2D Frame — Geometry + Loads + Reactions (Clean)');

xs = joints(:,2); ys = joints(:,3);
span = max(max(xs)-min(xs), max(ys)-min(ys));
if span == 0, span = 1; end

margin = 0.30 * span;
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

    plot([x1 x2], [y1 y2], '-', 'Color',[0.1 0.4 0.9], 'LineWidth',2);
end
plot(xs, ys, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);

% ===============================================================
% 4. PLOT SUPPORTS
% ===============================================================
if ~isempty(supports)
    for i = 1:size(supports,1)
        jid = supports(i,1);
        if ~isKey(id2row,jid), continue; end

        rr = id2row(jid);
        x  = joints(rr,2);
        y  = joints(rr,3);

        fixX = supports(i,2);
        fixY = supports(i,3);
        fixR = supports(i,4);

        if fixX==1 && fixY==1 && fixR==1
            drawFixedSupport(x,y,span);
        elseif fixX==1 && fixY==1 && fixR==0
            drawPinSupport(x,y,span);
        elseif (fixX==0 && fixY==1 && fixR==0) || (fixX==1 && fixY==0 && fixR==0)
            drawRollerSupport(x,y,span);
        else
            plot(x,y,'ks','MarkerFaceColor',[0.75 0.75 0.75],'MarkerSize',7);
        end
    end
end

% ===============================================================
% 5. PLOT JOINT LOADS (Fx/Fy separate + stacked "table"; Mz arc)
% ===============================================================
Lmax = 0.18 * span;
Lmin = 0.05 * span;

FxRef = 1; FyRef = 1; Mref = 1;

if ~isempty(jLoads)
    % Scaling refs for Fx/Fy separately (so each component scales nicely)
    FxAbs = abs(jLoads(:,2)); FxPos = FxAbs(FxAbs > tol);
    FyAbs = abs(jLoads(:,3)); FyPos = FyAbs(FyAbs > tol);

    if ~isempty(FxPos)
        FxRef = prctile(FxPos,75);
        if FxRef<=0, FxRef=max(FxPos); end
        if FxRef==0, FxRef=1; end
    end
    if ~isempty(FyPos)
        FyRef = prctile(FyPos,75);
        if FyRef<=0, FyRef=max(FyPos); end
        if FyRef==0, FyRef=1; end
    end

    % Moment scale reference
    Mabs = abs(jLoads(:,4)); Mpos = Mabs(Mabs > tol);
    if ~isempty(Mpos)
        Mref = prctile(Mpos,75);
        if Mref == 0, Mref = 1; end
    end

    for i = 1:size(jLoads,1)
        jid = jLoads(i,1);
        if ~isKey(id2row,jid)
            warning('J-Load references missing joint ID: %g', jid);
            continue;
        end

        rr = id2row(jid);
        x  = joints(rr,2);
        y  = joints(rr,3);

        Fx = jLoads(i,2);
        Fy = jLoads(i,3);
        Mz = jLoads(i,4);

        % ---- Fx arrow (global X), starting at the joint ----
        if abs(Fx) > tol
            Lx = Lmax * (abs(Fx)/(abs(Fx)+FxRef));
            Lx = max(Lx, Lmin);
            quiver(x, y, sign(Fx)*Lx, 0, 0, ...
                'Color',[0.95 0.75 0.1], 'LineWidth',2, 'MaxHeadSize',0.6);
        end

        % ---- Fy arrow (global Y), starting at the joint ----
        if abs(Fy) > tol
            Ly = Lmax * (abs(Fy)/(abs(Fy)+FyRef));
            Ly = max(Ly, Lmin);
            quiver(x, y, 0, sign(Fy)*Ly, 0, ...
                'Color',[0.95 0.75 0.1], 'LineWidth',2, 'MaxHeadSize',0.6);
        end

        % ---- Stacked joint force label (transparent table; no fill/border) ----
        % Each loaded joint gets its own table (your "if far away" requirement).
        if (abs(Fx) > tol) || (abs(Fy) > tol)
            txtF = sprintf('Fx = %.3g\nFy = %.3g', Fx, Fy);
            text(x, y - 0.12*span, txtF, ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', ...
                'FontSize',9, 'FontWeight','bold', ...
                'Color',[0.55 0.42 0.02], ...
                'BackgroundColor','none', 'EdgeColor','none', ...
                'Margin',1, 'Clipping','on');
        end

        % ---- Moment arc + label (kept) ----
        if abs(Mz) > tol
            r0 = 0.06*span + 0.06*span*(abs(Mz)/(abs(Mz)+Mref));
            drawMomentArc(x, y, r0, Mz, [0.85 0.2 0.55]);

            text(x + 1.10*r0, y + 0.20*r0, sprintf('Mz=%.3g',Mz), ...
                'Color',[0.55 0.05 0.30], 'FontSize',9, 'FontWeight','bold', ...
                'BackgroundColor','w', 'Margin',1, 'Clipping','on');
        end
    end
end

% ===============================================================
% 6. PLOT MEMBER LOAD MARKERS
% ===============================================================
if ~isempty(mLoads) && size(mLoads,2) >= 4
    markerLen = 0.06 * span;

    for k = 1:size(mLoads,1)
        mid  = mLoads(k,1);
        val  = mLoads(k,2);
        type = mLoads(k,3);
        ldir = mLoads(k,4);

        e = find(members(:,1)==mid,1);
        if isempty(e), continue; end

        r1 = id2row(members(e,2));
        r2 = id2row(members(e,3));

        x1 = joints(r1,2); y1 = joints(r1,3);
        x2 = joints(r2,2); y2 = joints(r2,3);

        Lm = hypot(x2-x1, y2-y1);
        if Lm < tol, continue; end

        if ldir == 1, gVec = [1;0]; else, gVec = [0;1]; end
        arrowSign = sign(val); if arrowSign == 0, arrowSign = 1; end

        if type == 1 % UDL markers
            dens = 7;
            for ii = 1:dens
                t = ii/(dens+1);
                xp = x1 + (x2-x1)*t;
                yp = y1 + (y2-y1)*t;

                quiver(xp, yp, gVec(1)*markerLen*arrowSign, gVec(2)*markerLen*arrowSign, 0, ...
                    'Color',[0.1 0.6 0.2], 'LineWidth',1.5, 'MaxHeadSize',0.7);
            end

            xm = (x1+x2)/2;
            ym = (y1+y2)/2;
            udlOff = 0.03*span;

            text(xm + gVec(1)*udlOff*arrowSign, ym + gVec(2)*udlOff*arrowSign, ...
                 sprintf('w=%.3g',val), 'Color',[0 0.45 0.15], ...
                 'FontSize',9, 'FontWeight','bold', ...
                 'BackgroundColor','w', 'Margin',1, 'Clipping','on');

        elseif type == 2 % midspan point load
            xp = (x1+x2)/2;
            yp = (y1+y2)/2;

            quiver(xp - gVec(1)*markerLen*arrowSign, yp - gVec(2)*markerLen*arrowSign, ...
                   gVec(1)*markerLen*arrowSign, gVec(2)*markerLen*arrowSign, 0, ...
                   'Color',[0 0.4 0.8], 'LineWidth',2.5, 'MaxHeadSize',0.8);

            text(xp + 0.01*span, yp + 0.02*span, sprintf('P=%.3g',val), ...
                'Color',[0 0.4 0.8], 'FontSize',9, 'FontWeight','bold', ...
                'BackgroundColor','w', 'Margin',1, 'Clipping','on');
        end
    end
end

% ===============================================================
% 7. PLOT SUPPORT REACTIONS (same as you had, but avoid forcing arrows if 0)
% ===============================================================
if nargin >= 6 && ~isempty(R) && ~isempty(supports)

    rxLen = 0.10*span;
    ryLen = 0.10*span;
    moRad = 0.07*span;

    for i = 1:size(supports,1)
        jid = supports(i,1);
        if ~isKey(id2row,jid), continue; end

        rr = id2row(jid);
        x  = joints(rr,2);
        y  = joints(rr,3);

        fixX = supports(i,2);
        fixY = supports(i,3);
        fixR = supports(i,4);

        fixedH = 0.065*span;
        [xa, ya] = supportCenterPoint(x,y,span,fixX,fixY,fixR,fixedH);

        Rx = R(3*rr-2);
        Ry = R(3*rr-1);
        Mz = R(3*rr);

        if fixX == 1 && abs(Rx) > tol
            quiver(xa - 0.06*span, ya + 0.01*span, sign(Rx)*rxLen, 0, 0, ...
                'Color',[0.85 0.1 0.1], 'LineWidth',2.2, 'MaxHeadSize',0.7);
        end

        if fixY == 1 && abs(Ry) > tol
            quiver(xa, ya - 0.10*span, 0, sign(Ry)*ryLen, 0, ...
                'Color',[0.85 0.1 0.1], 'LineWidth',2.2, 'MaxHeadSize',0.7);
        end

        if fixR == 1
            Mdraw = Mz;
            if abs(Mdraw) < tol, Mdraw = 1; end
            drawMomentArc(xa, ya, moRad, Mdraw, [0.85 0.1 0.1]);
        end

        RxDisp = Rx; if fixX~=1, RxDisp = 0; end
        RyDisp = Ry; if fixY~=1, RyDisp = 0; end
        MzDisp = Mz; if fixR~=1, MzDisp = 0; end

        text(xa, ya - 0.16*span, sprintf('Rx = %.3g\nRy = %.3g\nMz = %.3g', RxDisp, RyDisp, MzDisp), ...
            'HorizontalAlignment','center', 'VerticalAlignment','top', ...
            'FontSize',9, 'FontWeight','bold', ...
            'Color',[0.55 0 0], ...
            'BackgroundColor','none', 'EdgeColor','none', 'Margin',1, 'Clipping','on');
    end
end

hold off;
end