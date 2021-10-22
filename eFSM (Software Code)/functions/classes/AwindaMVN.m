classdef AwindaMVN
    properties
        Frames
        Frequency
        Data    % Struct Awinda.Data
        BodyWeight
        MSH_left
        MSH_right
        CoP
        Curve
    end
    methods
        %% determines the zero-lag filtered ground reaction force for all timesteps
        function F = getForce(obj)
            dt = 1/[obj.Frequency];  g = 9.81;
            %             v = diff(obj.Data.centerOfMass, [], 1)/dt;
            %             a = diff(v, [], 1)/dt;
            a = diff(obj.Data.centerOfMass, 2, 1)/dt^2;
            
            % F = - obj.BodyWeight * ( + a - g * repmat([0 0 1], size(a,1), 1));% ./ 9.81; %% m a = \sum {F_i} + m g
            F = obj.BodyWeight * ( + a + g * repmat([0 0 1], size(a,1), 1));% ./ 9.81; %% m a = \sum {F_i} + m g
            
            % butterworth zero lag filtering
            d = designfilt('lowpassiir','FilterOrder',1, ...
                'HalfPowerFrequency',0.09,'DesignMethod','butter'); % d1 = designfilt('lowpassfir', 'PassbandFrequency',0.09,'StopbandFrequency',0.5, 'PassbandRipple',1,'StopbandAttenuation',60, 'DesignMethod','equiripple');
            for i = 1:3
                %                 d = designfilt('lowpassiir','FilterOrder',12, ...
                %                     'HalfPowerFrequency',0.15,'DesignMethod','butter');
                F(:,i) = filtfilt(d,F(:,i));
            end
        end
        %% returns the angles of the total ground reaction force trajectory
        function [phi, psi, gamma] = getForceAngles(obj)
            Force = obj.getForce;
            phi = acos(Force(:,1)./vecnorm(Force')');
            psi = acos(Force(:,2)./vecnorm(Force')');
            gamma = acos(Force(:,3)./vecnorm(Force')');
        end
        %% returns calcaneus - , first metatarsal - and fifth metatarsal - landmarks of both feet
        function [LeftMarker, RightMarker] = getMarkerPoints(obj, i)
            % LFM, LVM, LCA
            LeftMarker = [obj.Data.position(i, 67:69) + obj.Data.segments.LeftToe.pLeftFirstMetatarsal([2 1 3])'; ...
                obj.Data.position(i, 67:69) + obj.Data.segments.LeftToe.pLeftFifthMetatarsal([2 1 3])'; ...
                obj.Data.position(i, 64:66) + obj.Data.segments.LeftFoot.pLeftHeelCenter([2 1 3])'];
            % RFM, RVM, RCA
            RightMarker = [obj.Data.position(i, 55:57) + obj.Data.segments.RightToe.pRightFirstMetatarsal([2 1 3])'; ...
                obj.Data.position(i, 55:57) + obj.Data.segments.RightToe.pRightFifthMetatarsal([2 1 3])'; ...
                obj.Data.position(i, 52:54) + obj.Data.segments.RightFoot.pRightHeelCenter([2 1 3])'];
            [RightMarker(1, :), RightMarker(2, :), ~] = Mirroring(RightMarker(2, :), RightMarker(1, :), RightMarker(3, :), ...
                obj.Data.position(i, 55:57) - [obj.Data.segments.RightToe.pRightToe([2 1])', -obj.Data.segments.RightToe.pRightToe(3)]);
        end
        %% draws the kinematics for the current state i
        function drawState(obj, i)
            State = reshape( obj.Data.position(i,:)', 3, length(obj.Data.position(i,:))/3 )';
            plot3(State(1:7,1), State(1:7,2), State(1:7,3),'r', ...        %pelvis -> spine -> head
                State(8:11,1), State(8:11,2), State(8:11,3),'r', ...       %right arm
                State(12:15,1), State(12:15,2), State(12:15,3),'r', ...    %left arm
                State(16:19,1), State(16:19,2), State(16:19,3),'r', ...    %right leg
                State(20:23,1), State(20:23,2), State(20:23,3),'r', ...    %left leg
                State(:,1),State(:,2),State(:,3),'b.', 'MarkerSize', 10,'HandleVisibility','off');
            hold on;
            plot3(State(7,1), State(7,2), State(7,3),'b.', 'MarkerSize', 50,'HandleVisibility','off');
            
            [LeftMarker, RightMarker] = obj.getMarkerPoints(i);
            % LFoot -> L HeelCenter -> LFM -> Toe -> LVM -> HeelCenter
            LFoot = [obj.Data.position(i, 64:66); ...
                LeftMarker(3,:); ...
                LeftMarker(1,:); ...
                obj.Data.position(i, 67:69);...
                LeftMarker(2,:); ...
                LeftMarker(3,:)];
            
            plot3(LFoot(:,1), LFoot(:,2), LFoot(:,3), 'r','HandleVisibility','off'); plot3(LFoot(:,1), LFoot(:,2), LFoot(:,3), 'k.','HandleVisibility','off');
            
            % RFoot -> R Heel -> RFM -> Toe -> RVM -> HeelCenter
            RFoot = [obj.Data.position(i, 52:54); ...
                RightMarker(3,:); ...
                RightMarker(1,:); ...
                obj.Data.position(i, 55:57); ...% + [sign(yawR)*obj.Data.segments.RightToe.pRightToe([2 1])', obj.Data.segments.RightToe.pRightToe(3)]; ...
                RightMarker(2,:); ...
                RightMarker(3,:)];
            plot3(RFoot(:,1), RFoot(:,2), RFoot(:,3), 'r','HandleVisibility','off'); plot3(RFoot(:,1), RFoot(:,2), RFoot(:,3), 'k.', 'MarkerSize', 10,'HandleVisibility','off');
            %             obj.drawPelvisOrientation(i);
        end
        %% Pelvis orientation
        % draws sagittal-axis in direction of movement
        function drawPelvisOrientation(obj, i)
            unit1 = [0,1,0,0];
            unit1 = quatTransform(unit1, obj.Data.orientation(i,1:4));
            unit1 = unit1(2:end);
            quiver3(obj.Data.position(i,1), obj.Data.position(i,2), 0*obj.Data.position(i,3), unit1(1), unit1(2), unit1(3));
        end
        % determines sagittal-axis, for time instance i, in direction of movement
        function [vec, q] = getPelvisOrientationVector(obj,i)
            q = quatTransform([0,1,0,0], obj.Data.orientation(i,1:4));
            vec = q(2:4);
        end
        % determines whole trajectory of sagittal-axis
        function vec = getTotalPelvisOrientationVector(obj)
            vec = zeros(obj.Frames, 3);
            for i = 1:obj.Frames
                vec(i, :) = obj.getPelvisOrientationVector(i);
            end
        end
        %% visualisation of the kinematic trajectory (a specified range can be declared)
        function showWholeAnimation(obj, range)
            %% play around
            F = obj.getForce;
            %%
            if nargin == 1
                range = 1:obj.Frames;
            end
            for i = range(1)+100 : range(end)-2%1:obj.Frames
                subplot(121)
                cla
                obj.drawState(i+2);
                axis equal
                %                 subplot(223)
                %                 cla
                %                 plot(F(i-100:i,1));
                %                 ylim([-200,200]); legend('F_{GRF,x}')
                %                 subplot(224)
                %                 cla
                %                 plot(F(i-100:i,2));
                %                 ylim([-200,200]); legend('F_{GRF,y}')
                subplot(122)
                cla
                plot(F(i-100:i,3));
                %                 ylim([-5000+obj.BodyWeight*9.81,5000+obj.BodyWeight*9.81]);  %legend('F_{GRF,z}')
                drawnow;
                %                 pause(0.01)
            end
        end
        %% function to compute the scalar nedded to project all segments onto the bottom plane
        function tr = get_to_bottom_transformed_points(obj, F_K)
            FootR = [52:54; 55:57]; FootL = [64:66; 67:69];
            LowArmL = [40:42; 43:45]; LowArmR = [28:30; 31:33];
            UpArmL = [37:39; 40:42]; UpArmR = [25:27; 28:30];
            LowLegL = [61:63; 64:66]; LowLegR = [49:51; 52:54];
            UpLegL = [58:60; 61:63]; UpLegR = [46:48; 49:51];
            
            t0 = 2;
            % % t = LineParam(F_K, point)
            % %             ->    t = point(:,3)./(point(:,3) - F_K(:,3));
            % % Feet
            tr.tFootL = ((obj.Data.position(1+t0:end,FootL(2,3)) + obj.Data.position(1+t0:end,FootL(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,FootL(2,3)) - obj.Data.position(1+t0:end,FootL(1,3)))./2) - F_K(:,3));
            tr.tFootR = ((obj.Data.position(1+t0:end,FootR(2,3)) + obj.Data.position(1+t0:end,FootR(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,FootR(2,3)) - obj.Data.position(1+t0:end,FootR(1,3)))./2) - F_K(:,3));
            % % Lower Arms
            tr.tLowArmL = ((obj.Data.position(1+t0:end,LowArmL(2,3)) + obj.Data.position(1+t0:end,LowArmL(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,LowArmL(2,3)) - obj.Data.position(1+t0:end,LowArmL(1,3)))./2) - F_K(:,3));
            tr.tLowArmR = ((obj.Data.position(1+t0:end,LowArmR(2,3)) + obj.Data.position(1+t0:end,LowArmR(1,3)))./2)./(-((obj.Data.position(1+t0:end,LowArmR(2,3)) - obj.Data.position(1+t0:end,LowArmR(1,3)))./2) - F_K(:,3));
            % % Upper Arms
            tr.tUpArmL = ((obj.Data.position(1+t0:end,UpArmL(2,3)) + obj.Data.position(1+t0:end,UpArmL(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,UpArmL(2,3)) - obj.Data.position(1+t0:end,UpArmL(1,3)))./2) - F_K(:,3));
            tr.tUpArmR = ((obj.Data.position(1+t0:end,UpArmR(2,3)) + obj.Data.position(1+t0:end,UpArmR(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,UpArmR(2,3)) - obj.Data.position(1+t0:end,UpArmR(1,3)))./2) - F_K(:,3));
            % % Upper Legs
            tr.tUpLegL = ((obj.Data.position(1+t0:end,UpLegL(2,3)) + obj.Data.position(1+t0:end,UpLegL(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,UpLegL(2,3)) - obj.Data.position(1+t0:end,UpLegL(1,3)))./2) - F_K(:,3));
            tr.tUpLegR = ((obj.Data.position(1+t0:end,UpLegR(2,3)) + obj.Data.position(1+t0:end,UpLegR(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,UpLegR(2,3)) - obj.Data.position(1+t0:end,UpLegR(1,3)))./2) - F_K(:,3));
            % % Shanks
            tr.tLowLegL = ((obj.Data.position(1+t0:end,LowLegL(2,3)) + obj.Data.position(1+t0:end,LowLegL(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,LowLegL(2,3)) - obj.Data.position(1+t0:end,LowLegL(1,3)))./2) - F_K(:,3));
            tr.tLowLegR = ((obj.Data.position(1+t0:end,LowLegR(2,3)) + obj.Data.position(1+t0:end,LowLegR(1,3)))./2)./(-0*((obj.Data.position(1+t0:end,LowLegR(2,3)) - obj.Data.position(1+t0:end,LowLegR(1,3)))./2) - F_K(:,3));
            % % Torso
            tr.tUpTorso = (obj.Data.position(1+t0:end,15))./(-0*obj.Data.position(1+t0:end,15) - F_K(:,3));
            tr.tMidTorso = (obj.Data.position(1+t0:end,9))./(-0*obj.Data.position(1+t0:end,9) - F_K(:,3));
            
            tr.tLowTorso = (obj.Data.position(1+t0:end,3))./(-0*obj.Data.position(1+t0:end,3) - F_K(:,3));
            tr.tHead = (obj.Data.position(1+t0:end,21))./(-0*obj.Data.position(1+t0:end,21) - F_K(:,3));
        end
        %% computes force shadows
        function [XY, xp, yp, f] = getForceShadow(obj, i, t0, F_K, tr, distribution, CoM, a, CoP, distanceFactor, Bool, extremities)
            minima = [min(obj.Data.position(i+t0,1:3:end)); min(obj.Data.position(i+t0,2:3:end)); min(obj.Data.position(i+t0,3:3:end))];
            maxima = [max(obj.Data.position(i+t0,1:3:end)); max(obj.Data.position(i+t0,2:3:end)); max(obj.Data.position(i+t0,3:3:end))];
            
            FootR = [52:54; 55:57]; FootL = [64:66; 67:69];
            LowArmL = [40:42; 43:45]; LowArmR = [28:30; 31:33];
            UpArmL = [37:39; 40:42]; UpArmR = [25:27; 28:30];
            LowLegL = [61:63; 64:66]; LowLegR = [49:51; 52:54];
            UpLegL = [58:60; 61:63]; UpLegR = [46:48; 49:51];
            LowTorso = 1:3; MidTorso = 7:9; UpTorso = 13:15; Head = 19:21;
            
%             xp = linspace(minima(1) - .35*abs(minima(1)), maxima(1) + .35*abs(maxima(1)), 50);
%             yp = linspace(minima(2) - .05*abs(minima(2)), maxima(2) + .05*abs(maxima(2)), 50);
            xp = linspace(obj.Data.position(i+t0,19)-0.5, obj.Data.position(i+t0,19)+.5, 50);
            yp = linspace(obj.Data.position(i+t0,20)-0.5, obj.Data.position(i+t0,20)+.5, 50);
            [X, Y] = meshgrid(xp,yp); XY = [X(:) Y(:)];
            
            RTorso = diag(a(1:2));
            if Bool.Use_Extremities
                RFootR = diag(extremities(1:2));
                RFootL = RFootR;
                RLowArmL = diag(extremities(3:4));
                RLowArmR = RLowArmL;
            end
            
            RUpArmL = diag(a(7-4:8-4)); RUpArmR = RUpArmL;
            RLowLegL = diag(a(9-4:10-4)); RLowLegR = RLowLegL;
            RUpLegL = diag(a(11-4:12-4)); RUpLegR = RUpLegL;
            CovTorso = obj.Two_Points_To_Covariance(obj.Data.position(i+t0,25:26)', obj.Data.position(i+t0,37:38)', RTorso);
            
            mu = obj.Data.position(i+t0, LowTorso(1:2)) - tr.tLowTorso(i)*(-F_K(i, 1:2) - 0*obj.Data.position(i+t0, LowTorso(1:2)));
            Maxs.MaxLowTorso = max(mvnpdf(XY, mu, CovTorso));
            % z = - (distribution.LowTorso/Maxs.MaxLowTorso)*mvnpdf(XY, mu, CovTorso);
            fLowTorso = @(x,y) - (distribution.LowTorso/Maxs.MaxLowTorso)*mvnpdf([x,y], mu, CovTorso);
            
            mu = obj.Data.position(i+t0, MidTorso(1:2)) - tr.tMidTorso(i)*(-F_K(i, 1:2) - 0*obj.Data.position(i+t0, MidTorso(1:2)));
            Maxs.MaxMidTorso = max(mvnpdf(XY, mu, CovTorso));
            % z = z - (distribution.MidTorso/Maxs.MaxMidTorso)*mvnpdf(XY, mu, CovTorso);
            fMidTorso = @(x,y) - (distribution.MidTorso/Maxs.MaxMidTorso)*mvnpdf([x,y], mu, CovTorso);
            
            mu = obj.Data.position(i+t0, UpTorso(1:2)) - tr.tUpTorso(i)*(-F_K(i, 1:2) - 0*obj.Data.position(i+t0, UpTorso(1:2)));
            Maxs.MaxUpTorso = max(mvnpdf(XY, mu, CovTorso));
            % z = z - (distribution.UpTorso/Maxs.MaxUpTorso)*mvnpdf(XY, mu, CovTorso);
            fUpTorso = @(x,y) - (distribution.UpTorso/Maxs.MaxUpTorso)*mvnpdf([x,y], mu, CovTorso);
            
            mu = obj.Data.position(i+t0, Head(1:2)) - tr.tHead(i)*(-F_K(i, 1:2) - 0*obj.Data.position(i+t0, Head(1:2)));
            Maxs.MaxHead = max(mvnpdf(XY, mu, eye(2)));
            % z = z - (distribution.Head/Maxs.MaxHead)*mvnpdf(XY, mu, eye(2));
            fHead = @(x,y) - (distribution.Head/Maxs.MaxHead)*mvnpdf([x,y], mu, eye(2));
            
%             for k = 5:length(string)
%                 range = eval(string(k));
%                 eval(strcat("R = R", string(k), ";"));
%                 eval(strcat("Max", string(k)," = max(mvnpdf(XY, .5.*(obj.Data.position(",num2str(i+t0),",", num2str(range(2,1)),":",num2str(range(2,2)),") + obj.Data.position(", ...
%                     num2str(i+t0),",", num2str(range(1,1)),":",num2str(range(1,2)),")) - tr.t", string(k), "(", num2str(i), ")*(-F_K(", num2str(i), ...
%                     ",1:2) - .5.*0*(obj.Data.position(",num2str(i+t0),",", num2str(range(2,1)),":",num2str(range(2,2)),") + obj.Data.position(",num2str(i+t0),",", num2str(range(1,1)), ...
%                     ":",num2str(range(1,2)),"))), obj.Two_Points_To_Covariance(obj.Data.position(", num2str(i+t0),",", num2str(range(1,1)), ...
%                     ":", num2str(range(1,3)), ")', obj.Data.position(", num2str(i+t0), ",", num2str(range(2,1)), ":", num2str(range(2,3)), ")', R)));"));
%                 
%                 eval(strcat("Maxs.Max", string(k), " = Max", string(k),";"));
%                 
%                 eval(strcat("f", string(k), " = @(x,y) -(distribution.", string(k), "/Max", string(k), ")*mvnpdf([x,y], .5.*(obj.Data.position(",num2str(i+t0),",", num2str(range(2,1)),":",num2str(range(2,2)),") + obj.Data.position(", ...
%                     num2str(i+t0),",", num2str(range(1,1)),":",num2str(range(1,2)),")) - tr.t", string(k), "(", num2str(i), ").*(-F_K(", num2str(i), ...
%                     ",1:2) - .5.*0*(obj.Data.position(",num2str(i+t0),",", num2str(range(2,1)),":",num2str(range(2,2)),") + obj.Data.position(",num2str(i+t0),",", num2str(range(1,1)), ...
%                     ":",num2str(range(1,2)),"))), obj.Two_Points_To_Covariance(obj.Data.position(", num2str(i+t0),",", num2str(range(1,1)), ...
%                     ":", num2str(range(1,3)), ")', obj.Data.position(", num2str(i+t0), ",", num2str(range(2,1)), ":", num2str(range(2,3)), ")', R));"));
%                 
%             end
            %%
            % k = 5
            range_temp = UpArmL;
            R = RUpArmL;
            MaxUpArmL = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tUpArmL(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxUpArmL = MaxUpArmL;
            fUpArmL = @(x,y) -(distribution.UpArmL/MaxUpArmL)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tUpArmL(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            % k = 6
            range_temp = UpArmR;
            R = RUpArmR;
            MaxUpArmR = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tUpArmR(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxUpArmR = MaxUpArmR;
            fUpArmR = @(x,y) -(distribution.UpArmR/MaxUpArmR)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tUpArmR(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            % k = 7
            range_temp = LowLegL;
            R = RLowLegL;
            MaxLowLegL = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tLowLegL(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxLowLegL = MaxLowLegL;
            fLowLegL = @(x,y) -(distribution.LowLegL/MaxLowLegL)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tLowLegL(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            % k = 8
            range_temp = LowLegR;
            R = RLowLegR;
            MaxLowLegR = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tLowLegR(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxLowLegR = MaxLowLegR;
            fLowLegR = @(x,y) -(distribution.LowLegR/MaxLowLegR)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tLowLegR(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            % k = 9
            range_temp = UpLegL;
            R = RUpLegL;
            MaxUpLegL = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tUpLegL(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxUpLegL = MaxUpLegL;
            fUpLegL = @(x,y) -(distribution.UpLegL/MaxUpLegL)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tUpLegL(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            % k = 10
            range_temp = UpLegR;
            R = RUpLegR;
            MaxUpLegR = max(mvnpdf(XY, .5.*(obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2))) ...
                - tr.tUpLegR(i)*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) + obj.Data.position(i+t0, range_temp(1,1) : range_temp(1,2)))), ...
                obj.Two_Points_To_Covariance(obj.Data.position(i+t0, range_temp(1,1) :range_temp(1,3))', obj.Data.position(i+t0, range_temp(2,1) : range_temp(2,3))', R)));
            Maxs.MaxUpLegR = MaxUpLegR;
            fUpLegR = @(x,y) -(distribution.UpLegR/MaxUpLegR)*mvnpdf([x,y], .5.*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2))) - tr.tUpLegR(i).*(-F_K(i,1:2) - .5.*0*(obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,2)) ...
                + obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,2)))), obj.Two_Points_To_Covariance(obj.Data.position(i+t0,range_temp(1,1) : range_temp(1,3))', obj.Data.position(i+t0,range_temp(2,1) : range_temp(2,3))', R));
            
            %%
            
            f = @(x,y) (fHead(x,y) + fLowTorso(x,y) + fMidTorso(x,y) + fUpTorso(x,y) + ...
                fUpArmL(x,y) + fUpArmR(x,y) + fLowLegL(x,y) + fLowLegR(x,y) + fUpLegL(x,y) + fUpLegR(x,y));
            
            if Bool.CoM_Correction
                f = @(x,y) (f(x + distanceFactor*(CoM(i+t0,1) - CoP(1)),y + distanceFactor*(CoM(i+t0,2) - CoP(2))) - f(maxima(1) + .25*abs(maxima(1)),maxima(2) + .25*abs(maxima(2))));  %% warum hab ich das gemacht?^^ -> Um "ausserhalb" Gewichtung 0 zu erhalten
            else
                f = @(x,y) (f(x,y) - f(maxima(1) + .25*abs(maxima(1)),maxima(2) + .25*abs(maxima(2))));  %% -> Um "ausserhalb" Gewichtung 0 zu erhalten
            end
        end
        function C = Two_Points_To_Covariance(obj, p1, p2, R)
            Diff = p1-p2;
            sgn=sign(p1(2)-p2(2));
            alpha = acos(Diff(1)/norm(p1-p2));
            C = [cos(alpha) -sgn*sin(alpha); sgn*sin(alpha) cos(alpha)]*R; %[1.6 0 ; 0 0.4]; %*[1.75 0 ; 0 0.25];
            C = C';
            
            C = (C'*C)./sqrt(norm(C));
        end
        %% heel strike tests for normal walking measurements
        function [ZeroIdxStrike, ZeroIdxOff] = HeelStrike(obj, side)
            %% Offline !
            switch side
                case "right"
                    HeelPos = obj.Data.position(:, 52:54) + obj.Data.segments.RightFoot.pRightHeelCenter([2 1 3])';
                    ToePos = obj.Data.position(:, 55:57);
                    Index = 54;
                case "left"
                    HeelPos = obj.Data.position(:, 64:66) + obj.Data.segments.LeftFoot.pLeftHeelCenter([2 1 3])';
                    ToePos = obj.Data.position(:, 67:69);
                    Index = 66;
            end
            PelvisPos = obj.Data.position(:,1:3);
            
            DistanceHeel = HeelPos - PelvisPos;
            DistanceToe = ToePos - PelvisPos;
            vec = obj.getTotalPelvisOrientationVector;
            differentialHeel = diff(sum(DistanceHeel.*vec,2),[],1)*60;
            differentialToe = diff(sum(DistanceToe.*vec,2),[],1)*60;
            
            ZeroIdxStrike = []; ZeroIdxOff = [];
            
            for i = 1:length(differentialHeel)-1
                if differentialHeel(i) > 0 && differentialHeel(i+1) < 0
                    ZeroIdxStrike(end+1) = i+1;
                end
                
                if differentialToe(i) < 0 && differentialToe(i+1) > 0
                    ZeroIdxOff(end+1) = i+1;
                end
            end
            
            % > 0.2 bei FL Load   nicht ben?tigt bei AlexW
            % 0.4 bei AV Load    (0.3 bei ENVISIBLE !)
            ZeroIdxStrike(find(sqrt(sum(obj.Data.velocity(ZeroIdxStrike, Index).^2,2)) > 0.3)) = [];
            %     I = find(obj.Data.position(ZeroIdxStrike, Index) > 0.2);
            %     ZeroIdxStrike(I) = [];
        end
    end
end

function [FM, VM, Toe] = Mirroring(FM, VM, Heel, Toe)
p_LCA = Heel;
n = cross(((FM+VM)/2 - Heel), [0 0 1]);
t = (p_LCA - VM)*n'/(n*n');
VM = VM + 2*t*n;

t = (p_LCA - FM)*n'/(n*n');
FM = FM + 2*t*n;

t = (p_LCA - Toe)*n'/(n*n');
Toe = Toe + 2*t*n;
end
