% the code calculates crevasse slpay morphology and water discharge outflow of the crevasse slpay
% that is triggered by a predesigned flood event and evolves afterwards.
% the codes are originally written for the Lower Yellow River(a suspended load dominated river) for the purpose of calculating sediment budget 
% on the floodplain over a long timescale,say as long as hundreds years, so the present module can not be applied to other alluvial rivers
% without modifying those lines related to channel geometry,bankfull discharge and bank erosion(depostion). 
% Please read the comments of the codes carefully and contact the author
% for further information on the model.Email: chenyunzhen1010@gmail.com

% inputs: Q-daily water discharge (m^3/s),1by1 matrix
%         S-sediment Cs (kg/m^3),1by1 matrix
%         delta_hcs0 (m) - the given initial depth of crevasse slpay ,usually set to 1m
%         Bcs0 (m) - the given initial width of crevasse slpay , usually set to 2m
%         Bc, delta_h, Bf, nc,nf,j - these are inputs for subroutine "flow.m"
%               Bc(m) - the width of the rectangular conceptual channel; 
%               delta_h(m) - the depth of the rectangular conceptual channel or the average depth in the cross section of channel;
%               (suppose the rectangular conceptual channel that has the same cross section area and full bank discharge as the real channel) )
%               Bf(m) - the width of floodplain(inside dikes);
%               nc - the manning coefficient of channel;
%               nf - the manning coefficient of floodplain;
%               j - the channel slope
%               
%         hs(m) - superelevation (the elevation of lowest point of channel bed, 
%              reference point for elevation: the elevation of the ground beyond the floodplain (no floodplain sedimentation )is set to 0)
%         Mb(kg/m^2/s) - M-coefficient for erosion rate for the bottom of crevasse slpay
%         Ms(kg/m^2/s) - M-coefficient for erosion rate for the two side slopes of crevasse slpay;  usually it is higher than Mb
%         ucr_e(m/s) - critical velocity for erosion 
%         ucr_d(m/s)- critical  velocity for depositon
%         dt(s) - time step (as Q is daily water discharge,the total calculating time is a day or 24*3600s) 
%         wdb - the width of dike at the root 
%         sg -  the slope of ground from reference point for elevation to the dike or the edge of floodplain
%         ws(m/s) - settling velocity of suspended load in the channel

% outputs: Delta_hcs - the depth of crevasse slpay at different times (time step dt) 
%          Bcs -  the width of crevasse slpay at different times(time step dt)
%          Zcs - the elevation of the bottom of CS at different times (time step dt)
%          Qcs - outflow discharge of CS at different time step at different times(time step dt)
%          Qout- water dischage to lower reach at different times (time step dt)
%          Hcs - water depth in the crevasse slpay at different times(time step dt)
%          Vcs - water velocity in the crevasse slpay at different times(time step dt)
%          Qd- averaged daily water dischage to lower reach of crevasse splay 
%          Delta_hcsd - averaged daily crevasse slpay depth
%          Bcsd - averaged daily crevasse slpay width
%          
% some variables in calculation:
%          t - the number of time steps druing a day
%          Qfb - bank full discharge
%          hmax - the biggest depth of the real channel
%          vci - flow velocity for Q in the rectangular conceptual channel;
%          hci -  averaged water depth for Q in the rectangular conceptual channel
%          vfi - flow velocity for Q on the rectangular conceptual floodplain;
%          hci - averaged water depth for Q on the floodplain;
%          wc - the width of flow for Q in the rectangular conceptual channel
%          wf - the width of flow for Q on the floodplain
%          Qci - water discharge in the rectangular conceptual channel;
%          hqc - the biggest water depth in the real channel
%          X,Ax -  two variables used in calculating hqc, see attached introduction file "introduction to modle Crevassesplay.doc" 
%          Ac - the cross section of flow in the rectangular conceptual channel
%          h0,ihqc,m,n,AR,AL,DA,dA,ii: the intermediate variables used in calculating hqc
%          Fi: Fraud nomber for the flow in the channel
%          Sv: volume sediment concentration 
%          L: the distance between reference point where elevation is 0 to the channel bank 
%          Zcsb: the bottom elevation of a crevasse slpay whose flow slope is equal to the channel slope j 
%          Qabove: the water discharge above the elevation of crevasseslpay
%                  bottom that can be distributed to the crevasseslpay
%                  outflow
%          rq: the discharge ratio of crevasseslpay outflow and Qabove
%          Xcs, Acs:  the intermediate variables used in calculating Qabove
%          jcs: the slope of the outflow of crevasse splay


function [Qd,Delta_hcs,Bcs,Delta_hcsd,Bcsd,Qout,Qcs,Zcs,Hcs,Vcs]=Crevassesplay(Q,S,delta_hcs0,Bcs0,Bc,delta_h,Bf,nc,nf,j,hs,Mb,Ms,ucr_e,ucr_d,dt,sg,wdb,ws)


% line 69-91 are for calculating the biggest water depth hqc for Q in the real channel
% The model uses a simplified conceptual method. If the real relationships between the biggest channel depth and water discharge are available, just use the real relationship   
Qfb=Bc*(delta_h)^(5/3)*sqrt(j)/nc;    % use Manning equation to calculate the bankfull discharge
hmax=sqrt(Bc)/34/Qfb^(-0.22);         % this is the relationship between the biggest depth of the real channel and the bankfull discharge 
                                      % specific for Yellow River (China).
                                      % Change this relationship if necessary
[vci,hci,vfi,hfi,wc,wf]=flow(Q,Bc,delta_h,Bf,nc,nf,j,1); % use Manning equation to calculate flow parameters in the corresponding rectangular conceptual channel  
Qci=vci*hci*wc;
X=delta_h*Bc/2/hmax;                  
Ax=X*(hmax-delta_h);
Ac=hci*wc;
if hci>delta_h                        
    hqc=hmax+(hci-delta_h);
elseif Ac<=delta_h*wc && Ac>Ax
    h0=hmax-delta_h;
    ihqc=h0:0.01:hmax;
    [m,n]=size(ihqc);
    AR=(X*ones(m,n)+(wc*delta_h*ones(m,n)-(hmax*ones(m,n)-ihqc)*(wc-2*X))/2/delta_h).*(ihqc-ones(m,n)*hmax+ones(m,n)*delta_h);
    AL=ones(m,n)*(Ac-Ax);
    DA=abs(AL-AR);
    [dA,ii]=min(DA);
    hqc=h0+(ii-1)*0.01;    
else hqc=sqrt(Ac*(hmax-delta_h)/X);      
end

t=round(24*3600/dt+1);      % the number of time steps druing a day
Delta_hcs=zeros(t+1,1); Bcs=zeros(t+1,1); Zcs=zeros(t+1,1); Qcs=zeros(t,1);Qout=zeros(t,1);Hcs=zeros(t,1);Vcs=zeros(t,1);
Delta_hcs(1,1)=delta_hcs0; Bcs(1,1)=Bcs0; Zcs(1,1)=hs+hmax-delta_hcs0;

Fi=Qci/sqrt(9.8*wc^2*hci^3);             % fraud no. (for calculating the discharge ratio of crevasse slpay outflow)

Sv=S/2650;                              

L=wdb+(hmax+hs)/sg;         % calculate the bottom elevation of the lowest point that a crevasse slpay is able to cut down   
Zcsb=L*j;         % rule1: the slope of crevasse slpay outflow is not less than the slope of lower channel
                  % rule2: the bottom elevation of a crevasse slpay is not lower than the the elevation of lowest point of channel bed(hs) 
                  % so the bottom elevation of the lowest point that a crevasse slpay is able to cut down  is max(hs,Zcsb)
for i=1:t
 % if the calculated Zcs is lower than max(hs,Zcsb), reset it to max(hs,Zcsb)    
if Zcs(i,1)<max(hs,Zcsb)             
   Zcs(i,1)=max(hs,Zcsb);
   Delta_hcs(i,1)=hmax+hs-Zcs(i,1);
end

% if the bottom elevation of a crevasse slpay is bigger than that of water suface in channel, there will be no outflow
if Zcs(i,1)>=hs+hqc         
    Qcs(i,1)=0;
    Qout(i,1)=Q;
    Hcs(i,1)=0;
    Vcs(i,1)=0;
    Zcs(i+1,1)=Zcs(i,1);
    Bcs(i+1,1)=Bcs(i,1);
    Delta_hcs(i+1,1)=Delta_hcs(i,1);
    
% if the bottom elevation of a crevasse slpay is smaller than that of water suface in channel and the crevasse slpay has not yet cut down to the lowest point    
% calculate the water dischage above the bottom of crevasse splay (Qabove) and the outflow of crevasse splay
% if the real relationship between Q and hqc is used, Qabove=Q-Q_bcs (Q_bcs is the water discharge when water surface elevation is equal to the bottom surface of crevasse splay)  
elseif Zcs(i,1)>max(hs,Zcsb) && Zcs(i,1)<hs+hqc     
                                                             
     if Zcs(i,1)<=hs+hmax-delta_h                                   
          Qabove=Q-vci*X*(Zcs(i,1)-hs)^2/(hmax-delta_h);             
      elseif Zcs(i,1)>hs+hmax-delta_h && Zcs(i,1)<=hs+hmax
          Xcs=(wc-(hs+hmax-Zcs(i,1))*(wc-2*X)/delta_h)/2;
          Acs=(X+Xcs)*(Zcs(i,1)-hs-hmax+delta_h);
          Qabove=Q-vci*(Acs+Ax);
      elseif Zcs(i,1)>hs+hmax
          Qabove=(wc*vci+wf*vfi)*(hs+hqc-Zcs(i,1));
      end
      
       if Qabove<0   
       Qabove=0;   
       end    
      
      rq=(1.55-1.45*Fi)*Bcs(i,1)/wc+0.16*(1-2*Fi);   % the discharge ratio of crevasse slpay outflow and Qabove, see "Wastewater Hydraulics:theory and practice"(Willi H Hager,1999,p505)
   if rq<=1
        Qcs(i,1)=rq*Qabove;
   else
        Qcs(i,1)=Qabove;
   end
   Qout(i,1)=Q-Qcs(i,1);    % water dischage to lower reach
   jcs=Zcs(i,1)/L;          % calculate the flow parameters for crevasse splay flow
   Hcs(i,1)=(nc*Qcs(i,1)/sqrt(jcs)/Bcs(i,1))^(3/5);
   Vcs(i,1)=Qcs(i,1)/Hcs(i,1)/Bcs(i,1);
   
   if Vcs(i,1)>ucr_e                                         % calculate the crevasse slpay morphology
           dEs=Ms*(Vcs(i,1)^2-ucr_e^2)/ucr_e^2*dt*2;         % now crevasse splay can both be widened and cut down
           Bcs(i+1,1)=Bcs(i,1)+dEs;
           dEb=Mb*(Vcs(i,1)^2-ucr_e^2)/ucr_e^2*dt;
           Zcs(i+1,1)=Zcs(i,1)-dEb;
           Delta_hcs(i+1,1)=Delta_hcs(i,1)+dEb;
    else
       if Vcs(i,1)<ucr_d && Vcs(i,1)>0
           dD=Sv*(1-Vcs(i,1)^2/ucr_d^2)*ws/0.6*dt;
           Zcs(i+1,1)=Zcs(i,1)+dD;
           Bcs(i+1,1)=Bcs(i,1);
           Delta_hcs(i+1,1)=Delta_hcs(i,1)-dD;
       else
           Zcs(i+1,1)=Zcs(i,1);
           Bcs(i+1,1)=Bcs(i,1);
           Delta_hcs(i+1,1)=Delta_hcs(i,1);
       end 
   end
% if the crevasse splay has cut down to the lowest point max(hs,Zcsb),crevasse splay can only be widened or silted vertically
elseif Zcs(i,1)==max(hs,Zcsb) && Zcs(i,1)<hs+hqc      
    if hs>=Zcsb                                       % if crevasse splay can cut down to the lowest point of channel bed, Qabove=Q
       rq=(1.55-1.45*Fi)*Bcs(i,1)/wc+0.16*(1-2*Fi);   
       if rq<=1
        Qcs(i,1)=rq*Q;
        else
        Qcs(i,1)=Q;
       end 
    elseif hs<Zcsb                              % if crevasse splay can not cut down to the lowest point of channel bed, calculate the water dischage above the bottom of crevasse splay and the outflow ration
        if Zcs(i,1)<=hs+hmax-delta_h                                  
          Qabove=Q-vci*X*(Zcs(i,1)-hs)^2/(hmax-delta_h);              
         elseif Zcs(i,1)>hs+hmax-delta_h && Zcs(i,1)<=hs+hmax         
          Xcs=(wc-(hs+hmax-Zcs(i,1))*(wc-2*X)/delta_h)/2;
          Acs=(X+Xcs)*(Zcs(i,1)-hs-hmax+delta_h);
          Qabove=Q-vci*(Acs+Ax);
          elseif Zcs(i,1)>hs+hmax
          Qabove=(wc*vci+wf*vfi)*(hs+hqc-Zcs(i,1));
        end
        
         if Qabove<0  
         Qabove=0;  
         end 
      
        rq=(1.55-1.45*Fi)*Bcs(i,1)/wc+0.16*(1-2*Fi);
       if rq<=1
        Qcs(i,1)=rq*Qabove;
       else
        Qcs(i,1)=Qabove;
       end
    end
   Qout(i,1)=Q-Qcs(i,1);
   jcs=Zcs(i,1)/L;
   Hcs(i,1)=(nc*Qcs(i,1)/sqrt(jcs)/Bcs(i,1))^(3/5);
   Vcs(i,1)=Qcs(i,1)/Hcs(i,1)/Bcs(i,1);
   
   if Vcs(i,1)>ucr_e                  % the crevasse splay has cut down to the lowest point max(hs,Zcsb),crevasse splay can only be widened or silted vertically
           dE=Ms*(Vcs(i,1)^2-ucr_e^2)/ucr_e^2*dt*2;
           Bcs(i+1,1)=Bcs(i,1)+dE;
           Zcs(i+1,1)=Zcs(i,1);
           Delta_hcs(i+1,1)=Delta_hcs(i,1);
   else
       if Vcs(i,1)<ucr_d && Vcs(i,1)>0
           dD=Sv*(1-Vcs(i,1)^2/ucr_d^2)*ws/0.6*dt;
           Zcs(i+1,1)=Zcs(i,1)+dD;
           Bcs(i+1,1)=Bcs(i,1);
           Delta_hcs(i+1,1)=Delta_hcs(i,1)-dD;
      
      else
       Zcs(i+1,1)=Zcs(i,1);
       Bcs(i+1,1)=Bcs(i,1);
       Delta_hcs(i+1,1)=Delta_hcs(i,1);
       
      end 
   end
   
end
 
end

Qd=mean(Qout);Delta_hcsd=mean(Delta_hcs(1:t,1));Bcsd=mean(Bcs(1:t,1)); % calculate the daily averaged water dischage to lower reach of crevasse splay and crevasse splay depth and width 



