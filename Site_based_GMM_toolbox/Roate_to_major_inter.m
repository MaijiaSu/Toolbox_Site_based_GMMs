function [acc_major,acc_inter] = Roate_to_major_inter(acc_H1,dt_H1,theta1,...
    acc_H2,dt_H2,theta2)
% Input:
% acc_H1 and acc_H2 are two horizontal time series
% dt_H1 and dt_H2 are sampling steps and should be eqaul
% theta1 and theta2 are azimuth angels (Colummns DL and DM in the NGA-west2 flatfile)

% if lengths of acc_H1 and acc_H2 are different, zero pad the shortest
% before rotating
if length(acc_H1)>length(acc_H2)
    acc_H2=[acc_H2; zeros(length(acc_H1)-length(acc_H2),1)];
elseif length(acc_H2)>length(acc_H1)
    acc_H1=[acc_H1; zeros(length(acc_H2)-length(acc_H1),1)];
end

% Compute the angle different
diff = theta1-theta2+360;
diff = rem(diff,360);
diff(diff==270) = -90;
if diff == -90
    tempacc_H1 = acc_H2;
    tempacc_H2 = acc_H1;
    temptheta1 = theta2;
    temptheta2 = theta1;
else %    if diff(i)==90 keep the same, also for i=34 and i=45
    tempacc_H1 = acc_H1;
    tempacc_H2 = acc_H2;
    temptheta1 = theta1;
    temptheta2 = theta2;
end

%% Preliminary calculations
if dt_H1 == dt_H2
    v1=sum(tempacc_H1.^2)*dt_H1;
    v2=sum(tempacc_H2.^2)*dt_H2;
    rho1_2=corr(tempacc_H1,tempacc_H2);
else
    error('dt_H1~=dt_H2');
end

%% Rotation angle (ccw) from Acc 1 and Acc 2 to major and intermediate
theta_rot=0.5*atan((2*rho1_2.*sqrt(v1).*sqrt(v2))./(v1-v2));
theta_rot(theta_rot<0)=theta_rot(theta_rot<0)+pi/2;
theta_rot_deg=theta_rot/pi*180;

%% acc_H1 and acc_H2 rotated to major and intermediate (major has largest AI)
% angle between major and Inter is either +90 or -90 degrees
Acc_rot1=cosd(theta_rot_deg)*tempacc_H1+sind(theta_rot_deg)*tempacc_H2;
Acc_rot2=-sind(theta_rot_deg)*tempacc_H1+cosd(theta_rot_deg)*tempacc_H2;

AI_rot1=sum(Acc_rot1.^2);
AI_rot2=sum(Acc_rot2.^2);

if AI_rot1>AI_rot2
    acc_major=Acc_rot1;
    acc_inter=Acc_rot2;
else
    acc_major=Acc_rot2;
    acc_inter=Acc_rot1;
end

end