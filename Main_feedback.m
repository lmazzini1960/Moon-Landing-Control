%% determine the extremal of the moon landing problem
close all
clear all
global x0 A0 g vref
% x0(1) initial altitude  m
% x0(2) initial  horizontal velocity  m/s
% x0(3) initial  vertical velocity     m/s
% tf   final time
% t time at which we output the solution

%% thrust
g=1.62; % moon gravitational accelleration in Pontani
A0=6*g; % start accelleraion;

%% load feedback
% load Landing_K_2D_C_MEM3g.mat vxr vyr alfa0 k0 V Crash
load Landing_K_2D_C_MEM6g.mat vxr vyr alfa0 k0 V Crash

%% extremal state solution
x0=zeros(3,1);
q0=zeros(3,1);
x0(1)=15000;
x0(2)=1692; 
x0(3)=0;  
a0=180;      % guess of initial firing angle
af=150;      % guess of the final firing angle
q0(1)=(-A0/x0(2))*(asinh(tand(af))-asinh(tand(a0)));
q0(2)=tand(a0);
q0(3)=(-x0(2)/A0)* (tand(a0)-tand(af))/( asinh(tand(af))-asinh(tand(a0)) );
options = optimset('TolX',1e-10,'TolFun',1e-6);
vref=0.01*x0(2);
[q,fval]=fminsearch(@finalcost,q0,options);
a0=atand(q(2))+180;
af=atand( q(2)-q(1)*q(3) )+180;
tf_ext=q(3);
%% calculate state evolution of the extremal
landing=0;
Nt_ext=2001;
t_ext=linspace(0,tf_ext,Nt_ext)';
x_ext=zeros(3,Nt_ext);
vxr_ext=zeros(Nt_ext,1);
vyr_ext=zeros(Nt_ext,1);
alfa_ext=zeros(Nt_ext,1);
for i=1:Nt_ext
x_ext(1,i)=x0(1)+( x0(3)- (A0/q(1))*sqrt(1+q(2)^2) )*t_ext(i)-g*t_ext(i)^2/2....
  -(A0/q(1)^2)*( (q(2)-q(1)*t_ext(i))/2 * sqrt(1+(q(2)-q(1)*t_ext(i))^2) )....
  -(A0/(2*q(1)^2))* asinh(q(2)-q(1)*t_ext(i)) + (A0/q(1)^2)*( (q(2)/2)*sqrt(1+q(2)^2)+0.5*asinh(q(2)) );  
x_ext(2,i)=x0(2)+(A0/q(1))*(asinh(q(2)-q(1)*t_ext(i))-asinh(q(2)));
x_ext(3,i)=x0(3)-g*t_ext(i)+(A0/q(1))*( sqrt(1+(q(2)-q(1)*t_ext(i))^2)- sqrt(1+q(2)^2) );
alfa_ext(i)=atand(q(2)-q(1)*t_ext(i))+180;
vxr_ext(i,1)=x_ext(2,i)/sqrt(A0*x_ext(1,i));
vyr_ext(i,1)=x_ext(3,i)/sqrt(A0*x_ext(1,i));
if x_ext(1,i) < 0.5
        landing=1;
end
    if landing==1
       vxr_ext(i,1)= vxr_ext(i-1,1);
       vyr_ext(i,1)= vyr_ext(i-1,1);
    end
end
final_error=sqrt( (x_ext(1,Nt_ext)/x0(1))^2+ (x_ext(2,Nt_ext)/vref)^2+ (x_ext(3,Nt_ext)/vref)^2 )
result=q'

%% determine tfinal from feedback table
vxr0=x0(2)/sqrt(A0*x0(1));
vyr0=x0(3)/sqrt(A0*x0(1));
[vx,ivx]=min(abs(vxr0-vxr));
[vy,ivy]=min(abs(vyr0-vyr));
%% determine landing duration
tfinal_feed=( V(ivy,ivx)+ (vxr0-vxr(ivx,1))*(V(ivy,ivx+1)-V(ivy,ivx))/(vxr(ivx+1,1)-vxr(ivx,1)) +...
        (vyr0-vyr(ivy,1))*(V(ivy+1,ivx)-V(ivy,ivx))/(vyr(ivy+1,1)-vyr(ivy,1)) ) * sqrt(x0(1)/A0);
%% initialize feedback vectors
% Nt_feed=3001;
% t_feed=linspace(0,tfinal_feed,Nt_feed)';
% r_feed=zeros(Nt_feed,1);
% vx_feed=zeros(Nt_feed,1);
% vy_feed=zeros(Nt_feed,1);
% vxr_feed=zeros(Nt_feed,1);
% vyr_feed=zeros(Nt_feed,1);
% x_feed=zeros(3,Nt_feed);
% alfa_feed=zeros(Nt_feed,1);

tau=tfinal_feed/5000;

%% initialize state
r_feed(1,1)=x0(1);
vx_feed(1,1)=x0(2);
vy_feed(1,1)=x0(3);
vxr_feed(1,1)=vxr0;
vyr_feed(1,1)=vyr0;
t_feed(1,1)=0;
tf_feed(1,1)=tfinal_feed;
%% simulation flags
landing=0;        % landing flag
pre_landing=0;    % pre-landing phase
failure=0;        % failure  phase 
r_landing=500;
r_start_landing=1000;
r_start_failure=8000; % 
r_end_failure=7000;
time_failure=0; % reset failure time counter
failure_duration=10; % duration of failure=recovery time
failure_angle=-45*pi/180;% max rotation in failure duration
% r_start_failure=0;
% r_end_failure=0;
%% loop with feedback table
i=2;
 while tf_feed(i-1,1) > 0.1   % land with tf_feed =0 alternatively use r_feed(i-1,1) > 1
     
    if r_feed(i-1,1) < r_start_failure
        failure=1;
    end
    % if r_feed(i-1,1) < r_end_failure   % for stop firing failure
    if time_failure >2* failure_duration  % for rotation failure
        failure=0;
    end
        
    if r_feed(i-1,1) < r_start_landing
        landing=0;
        pre_landing=1;
    end
    
    if r_feed(i-1,1) < r_landing
        landing=1;
        pre_landing=0;
    end
    
    if pre_landing==1  % in this phase we memorize the variation of the vxr,vyr but use the vxr,vyr table to feedback
        vxr_feed(i-1,1)=vx_feed(i-1,1)/sqrt(A0*r_feed(i-1,1));
        vyr_feed(i-1,1)=vy_feed(i-1,1)/sqrt(A0*r_feed(i-1,1));
        dvxr=(vxr_feed(i-1,1)- vxr_feed(i-2,1));
        dvyr=(vyr_feed(i-1,1)- vyr_feed(i-2,1));
    end

    if landing==1   % in this phase we use the variations of vxr,vyr  to feedback
       vxr_feed(i-1,1)= vxr_feed(i-2,1)+dvxr;
       vyr_feed(i-1,1)= vyr_feed(i-2,1)+dvyr;
    end
    
    if landing==0   % in this phase we  use the vxr,vyr table to feedback
        vxr_feed(i-1,1)=vx_feed(i-1,1)/sqrt(A0*r_feed(i-1,1));
        vyr_feed(i-1,1)=vy_feed(i-1,1)/sqrt(A0*r_feed(i-1,1));   
    end
    
    [gevxr(i-1,1),ivx]=min(abs(vxr_feed(i-1,1)-vxr));
    [gevyr(i-1,1),ivy]=min(abs(vyr_feed(i-1,1)-vyr));
    
%% determine landing duration feedback ad duration

    if ivx<length(vxr)
    alfa_feed(i-1,1) = alfa0(ivy,ivx) +  (vxr_feed(i-1,1)-vxr(ivx,1))*(alfa0(ivy,ivx+1)-alfa0(ivy,ivx))/(vxr(ivx+1,1)-vxr(ivx,1)); 
    else
    alfa_feed(i-1,1) = alfa0(ivy,ivx) + ( vxr_feed(i-1,1)-vxr(ivx,1) )*( alfa0(ivy,ivx)-alfa0(ivy,ivx-1) )/( vxr(ivx,1)-vxr(ivx-1,1) ); 
    end
    if ivy<length(vyr)
    alfa_feed(i-1,1) = alfa_feed(i-1,1)+ (vyr_feed(i-1,1)-vyr(ivy,1))*(alfa0(ivy+1,ivx)-alfa0(ivy,ivx))/(vyr(ivy+1,1)-vyr(ivy,1));  
    else
    alfa_feed(i-1,1) = alfa_feed(i-1,1)+ (vyr_feed(i-1,1)-vyr(ivy,1))*(alfa0(ivy,ivx)-alfa0(ivy-1,ivx))/(vyr(ivy,1)-vyr(ivy-1,1));  
    end

    if landing==0
    if ivx<length(vxr)
    tf_feed(i,1)   = V(ivy,ivx)*sqrt(r_feed(i-1)/A0) - tau + sqrt(r_feed(i-1)/A0) * (vxr_feed(i-1,1)-vxr(ivx,1))*(V(ivy,ivx+1)-V(ivy,ivx))/(vxr(ivx+1,1)-vxr(ivx,1)); 
    else
    tf_feed(i,1)   = V(ivy,ivx)*sqrt(r_feed(i-1)/A0) - tau + sqrt(r_feed(i-1)/A0) * ( vxr_feed(i-1,1)-vxr(ivx,1) )*( V(ivy,ivx)-V(ivy,ivx-1) )/( vxr(ivx,1)-vxr(ivx-1,1) ) ;
    end
    if ivy<length(vyr)
    tf_feed(i,1) = tf_feed(i,1) + sqrt(r_feed(i-1)/A0) * (vyr_feed(i-1,1)-vyr(ivy,1))*(V(ivy+1,ivx)-V(ivy,ivx))/(vyr(ivy+1,1)-vyr(ivy,1));
    else
    tf_feed(i,1) = tf_feed(i,1)  + sqrt(r_feed(i-1)/A0) * (vyr_feed(i-1,1)-vyr(ivy,1))*(V(ivy,ivx)-V(ivy-1,ivx))/(vyr(ivy,1)-vyr(ivy-1,1));
    end
    else % during landing phase final time is fixed
    tf_feed(i,1) = tf_feed(i-1,1) - tau;
    end

    %% dynamics in Euler format
    if failure==0  % in this phase we  feedback
    t_feed(i,1)=t_feed(i-1,1)+tau;
    r_feed(i,1)=r_feed(i-1)+vy_feed(i-1,1)*tau;
    vx_feed(i,1)=vx_feed(i-1)+A0*cos(alfa_feed(i-1,1))*tau;
    vy_feed(i,1)=vy_feed(i-1)-g*tau + A0*sin(alfa_feed(i-1,1))*tau;
    else           % in this phase we do not feedback
    time_failure=time_failure+tau;
    if time_failure < failure_duration
    rotation=failure_angle*(time_failure/failure_duration); % 30deg max rotation
    else
    rotation=failure_angle*(1-(time_failure-failure_duration)/failure_duration) ;
    end
    t_feed(i,1)=t_feed(i-1,1)+tau;
    r_feed(i,1)=r_feed(i-1)+vy_feed(i-1,1)*tau;
    vx_feed(i,1)=vx_feed(i-1)+A0*cos(alfa_feed(i-1,1)+rotation)*tau;
    vy_feed(i,1)=vy_feed(i-1)-g*tau+ A0*sin(alfa_feed(i-1,1)+rotation)*tau;   
    end
    %% update arrays counter
     i=i+1;
end
vxr_feed(i-1,1)= vxr_feed(i-2,1);
vyr_feed(i-1,1)= vyr_feed(i-2,1);
alfa_feed(i-1,1)=alfa0(ivy,ivx);
gevxr(i-1,1)=gevxr(i-2,1);
gevyr(i-1,1)=gevyr(i-2,1);

%% set plot variables
x_feed(1,:)=r_feed; x_feed(2,:)=vx_feed; x_feed(3,:)=vy_feed;
alfa_feed=alfa_feed*180/pi;

%% plot solution
figure(2)
plot(t_ext,x_ext(1,:),'k--',t_feed,x_feed(1,:),'r.','LineWidth',3)
title('altitude  [m]')
xlabel('[sec]')
ylabel('[m]')
grid on

figure(3)
plot(t_ext,x_ext(2,:),'k--',t_feed,x_feed(2,:),'r.','LineWidth',3)
title('horizontal speed [m/s]')
xlabel('[sec]')
ylabel('[m/s]')
grid on

figure(4)
plot(t_ext,x_ext(3,:),'k--',t_feed,x_feed(3,:),'r.','LineWidth',3)
title('vertical speed [m/s]')
xlabel('[sec]')
ylabel('[m/s]')
grid on

figure(5)
plot(t_ext,vxr_ext,'k--',t_feed,vxr_feed,'r.','LineWidth',3)
title('vx/sqrt(A0 r)')
xlabel('[sec]')
ylabel('[]')
grid on

figure(6)
plot(t_ext,vyr_ext,'k--',t_feed,vyr_feed,'r.','LineWidth',3)
title('vy/sqrt(A0 r)')
xlabel('[sec]')
ylabel('[]')
grid on

figure(7)
plot(t_ext,alfa_ext,'k--',t_feed,alfa_feed,'r.','LineWidth',3)
title('alfa')
xlabel('[sec]')
ylabel('[deg]')
grid on

figure(8)
plot(t_feed,tf_feed,'r.','LineWidth',2)
title('current estimated final time tf')
xlabel('[sec]')
ylabel('[sec]')
grid on

figure(9)
plot(t_feed,gevxr(:,1),'k',t_feed,gevyr(:,1),'r.','LineWidth',2)
title('grid error vxr/vyr')
grid on

R=1700000.0; %Moon Radius
ares=zeros(2,Nt_ext);
for i=1:Nt_ext
 ares(1,i)=x_ext(2,i)*x_ext(3,i)/(R+x_ext(1,i));
 ares(2,i)=x_ext(2,i)^2/(R+x_ext(1,i));
end    
  
figure(10)
plot(t_ext,ares(1,:),'k--',t_ext,ares(2,:),'r.','LineWidth',3)
title('non linear accellerations [m/s^2]')
xlabel('[sec]')
ylabel('[m/s^2]')
grid on

return


function fun=finalcost(q)
global x0 A0 g vref
%% calculates the function to be minimized to provide the correct costate
%% q are the minimization parameters 
% x(1)=y          x(2)=vx       x(3)=vy
% q(1)=lambda1'  q(2)=lambda4'  q(3)=tf
yf=x0(1)+( x0(3)- (A0/q(1))*sqrt(1+q(2)^2) )*q(3)-g*q(3)^2/2....
  -(A0/q(1)^2)*( (q(2)-q(1)*q(3))/2 * sqrt(1+(q(2)-q(1)*q(3))^2) )....
  -(A0/(2*q(1)^2))* asinh(q(2)-q(1)*q(3)) + (A0/q(1)^2)*( (q(2)/2)*sqrt(1+q(2)^2)+0.5*asinh(q(2)) );  
vxf=x0(2)+(A0/q(1))*(asinh(q(2)-q(1)*q(3))-asinh(q(2)));
vyf=x0(3)-g*q(3)+(A0/q(1))*( sqrt(1+(q(2)-q(1)*q(3))^2)- sqrt(1+q(2)^2) );
fun=sqrt((yf/1000)^2+(vxf/vref)^2+(vyf/vref)^2 );
end
