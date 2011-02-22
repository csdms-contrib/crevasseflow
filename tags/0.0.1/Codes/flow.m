% this is a code for the calculation of flow parameters using Manning equation, suppose a rectangular conceptual channel 
% inputs:
%   it- the number of water discharge inputs 
%   Q(m^3/s)-water dischage, it by 1 matrix
%   Bc(m)-the channel width when water discharge is the bank full discharge
%   Bf(m)-the width of floodplain(within dikes)
%   delta_h(m)-average water depth when water discharge is the bankfull discharge, it by 1 matrix, and it can be calculated with manning equation
%   nc-manning coefficient of channel
%   nf-manning coefficient of floodplain
%   j-channel slope
% outputs:
%   vc(m/s)- flow velocity in the channel 
%   vf(m/s)-flow velocity on the floodplain
%   wc(m)-water surface width in the channel
%   wf(m)- water surface width on the floodplain



function [vc,hc,vf,hf,wc,wf]=flow(Q,Bc,delta_h,Bf,nc,nf,j,it)
hc=(nc*Q/sqrt(j)/Bc).^(3/5);
hf=zeros(it,1);
wc=Bc*ones(it,1);
wf=zeros(it,1);
for t=1:it
   if hc(t,1)>delta_h(t,1)       % the flow spills over the floodplain
   ih=0.01:0.01:3;
   Qcal=1/nf*sqrt(j)*(Bf*ih).^(5/3)./(Bf*ones(1,300)+ih).^(2/3)+1/nc*sqrt(j)*(Bc*(delta_h(t,1)*ones(1,300)+ih)).^(5/3)./((Bc+2*delta_h(t,1))*ones(1,300)+ih).^(2/3);
   DQ=abs(Q(t,1)*ones(1,300)-Qcal);
   [dQ,i]=min(DQ);
   hf(t,1)=i*0.01;
   hc(t,1)=delta_h(t,1)+hf(t,1);
   wf(t,1)=Bf;
   end
end
vc=1/nc*hc.^(1/6).*(hc*j).^0.5;
vf=1/nf*hf.^(1/6).*(hf*j).^0.5;

    

