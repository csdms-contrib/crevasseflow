% this is a sample code for  crevasse splay modelling
%inputs: 
%    it- the number of input Q  
%    Q(m^3/s) - daily water discharge series, it by 1 matrix
%    Qs(kg/s) - daily sediment flux series, it by 1 matrix
%    ics- the time when crevasse splay begins
%    delta_hcs0,Bcs0,Bc,delta_h,Bf,nc,nf,j,hs,Mb,Ms,ucr_e,ucr_d,dt,sg,wdb,ws
%       -the inputs of Crevassesplay.m
%       for example, for modern Yellow River:
%       delta_hcs0=1,Bcs0=2,Bc=965,delta_h=2.3,Bf=4795,nc=0.009,nf=0.03,j=1.377e-4,hs=3,Mb=0.0005,Ms=0.004,ucr_e=1.5,ucr_d=0.7,dt=1800,sg=2.5e-4,wdb=25,ws=5e-4
% outputs:
%     Qd (m/s) -averaged daily water dischage to lower reach of crevasse splay, it by 1 matrix
%     Delta_hcs(m) - the depth of crevasse slpay at different times (time step dt),(24*3600/dt+2) by it by 1 matrix
%     Bcs(m) -  the width of crevasse slpay at different times(time step dt),(24*3600/dt+2) by it by 1 matrix
%     Zcs - the elevation of the bottom of CS at different times (time step dt), (24*3600/dt+2) by it by 1 matrix
%     Qcs - outflow of CS at different times(time step dt), (24*3600/dt+1) by it by 1 matrix
%     Qout- water dischage to lower reach at different times (time step dt), (24*3600/dt+1) by it by 1 matrix
%     Hcs - water depth in the crevasse slpay at different times(time step dt), (24*3600/dt+1) by it by 1 matrix
%     Vcs - water velocity in the crevasse slpay at different times(time step dt), (24*3600/dt+1) by it by 1 matrix
%     Delta_hcsd - averaged daily crevasse slpay depth, it by 1 matrix
%     Bcsd - averaged daily crevasse slpay width, it by 1 matrix
% As a test, you can use the attached input data Q and Qs, and run 
% [Qd,Delta_hcs,Bcs,Delta_hcsd,Bcsd,Qout,Qcs,Zcs,Hcs,Vcs]=mainCS(Q,Qs,365,232,1,2,965,2.3,4795,0.009,0.03,1.377e-4,-2,0.0005,0.004,1.5,0.7,1800,2.5e-4,25,4.5e-4);

function [Qd,Delta_hcs,Bcs,Delta_hcsd,Bcsd,Qout,Qcs,Zcs,Hcs,Vcs]=mainCS(Q,Qs,it,ics,delta_hcs0,Bcs0,Bc,delta_h,Bf,nc,nf,j,hs,Mb,Ms,ucr_e,ucr_d,dt,sg,wdb,ws)
S=Qs./Q;
t=round(24*3600/dt+1);
Qd=zeros(it,1);Delta_hcsd=zeros(it,1);Bcsd=zeros(it,1);Qout=zeros(t,1,it);Qcs=zeros(t,1,it);Delta_hcs=zeros(t+1,1,it);Bcs=zeros(t+1,1,it);
Zcs=zeros(t+1,1,it);Hcs=zeros(t,1,it);Vcs=zeros(t,1,it);
for i=1:ics-1
   Qd(i,1)=Q(i,1);Qout(:,1,i)=Q(i,1);
end
[Qd(ics,1),Delta_hcs(:,1,ics),Bcs(:,1,ics),Delta_hcsd(ics,1),Bcsd(ics,1),Qout(:,1,ics),Qcs(:,1,ics),Zcs(:,1,ics),Hcs(:,1,ics),Vcs(:,1,ics)]=Crevassesplay(Q(ics,1),S(ics,1),delta_hcs0,Bcs0,Bc,delta_h,Bf,nc,nf,j,hs,Mb,Ms,ucr_e,ucr_d,dt,sg,wdb,ws);
for i=ics+1:it
    [Qd(i,1),Delta_hcs(:,1,i),Bcs(:,1,i),Delta_hcsd(i,1),Bcsd(i,1),Qout(:,1,i),Qcs(:,1,i),Zcs(:,1,i),Hcs(:,1,i),Vcs(:,1,i)]=Crevassesplay(Q(i,1),S(i,1),Delta_hcs(t+1,1,i-1),Bcs(t+1,1,i-1),Bc,delta_h,Bf,nc,nf,j,hs,Mb,Ms,ucr_e,ucr_d,dt,sg,wdb,ws);
end
T=1:it;
plot(T,Q,'-k',T,Qd,'-r')


    
