clc
clear all
close all
L=100;
preamble=[1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1; ...
1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;-1; ...
1;-1;1;-1;1;-1;1;-1;1;-1;1;-1; ...
1;-1;1;-1;1;-1;-1;1;-1;1;-1;1;-1;1;1;-1;1;-1;-1;1;-1;1;1;-1];
snr=10;
sims=500000;
nbits=120;
nprea=88;
pithsymbol=zeros(length(snr),1);
for i=1:length(snr)
    biterror=0;
    for j=1:sims        
        [bits,y_miller] = RN_miller(nbits);    %dimiourgia 20 bits kai ta kolame meta to preamble
        x = [preamble; y_miller];
        x_upsample = upsample(x,L/4);
        %tora to t*
        h= rayleigh_fading(true);
        metro=abs(h)^2;
        sigma = sqrt( metro*L/(2*10^(snr(i)/10)));
        y_conv = conv(x_upsample(1:end),ones(L/4,1));
        n = sigma*randn(size(y_conv)) + 1i*sigma*randn(size(y_conv));
        y_r = h * y_conv + n;
        y_r_filt=conv(y_r,ones(L/4,1));
        y_r_sampled  = y_r_filt(L/4:L/4:length(y_conv));
        y_info=y_r_sampled(length(preamble)+1:end);
        %tora to t*
%         for iii=0:L
%             sum=0;
%             for jj=1:nprea*L/4
%                 sum=sum+x_upsample(jj)*y_r_filt(jj+iii);
%             end
%             t(iii+1)=abs(sum);
%         end
%         [~,tau] = max(t);
%         %tora ektimisi tou h
%         sum=0;
%         for iii=1:nprea*L/4
%             sum=sum+x_upsample(iii)^2;
%         end
%         sum1=0;
%         for ii=tau+1:nprea*L/4+tau
%             sum1=sum1+x_upsample(ii-tau)*y_r_filt(ii);
%         end
%         hnew=sum1/sum;
%         h_e=4*hnew/L;
        h_e=h;
        %tora tha kanoume detect viterbi
        detect=zeros(nbits,1);
        p=.25; %arxiki katanomi ton simvolon
        dmax=zeros(4,1);
        A=1000;
        dd=zeros(4,nbits);
        sigma2=sigma*sqrt(L/2);
        for ii=1:4:4*nbits-3
            y0=y_info(ii);
            y1=y_info(ii+1);
            y2=y_info(ii+2);
            y3=y_info(ii+3);
            %
            %2 ta prosarmosmena filtra
            r1=(y0-y1+y2-y3)/(2);
            r2=(y0-y1-y2+y3)/(2);
            %
            if(ii==1)
             dmax(1)=log(p*gauss(r1-2*h_e,r2,sigma2)); %S1
             dmax(2)=log(p*gauss(r1+2*h_e,r2,sigma2)); %s2
             dmax(3)=log(p*gauss(r2+2*h_e,r1,sigma2)); 
             dmax(4)=log(p*gauss(r2-2*h_e,r1,sigma2)); 
             tempd1=dmax(1);
             tempd2=dmax(2);
             tempd3=dmax(3);
             tempd4=dmax(4);
             dd(1,1)=1;
             dd(2,1)=2;
             dd(3,1)=3;
             dd(4,1)=4;
            else   
             temp1=log(gauss(r1-2*h_e,r2,sigma2));
             temp2=log(gauss(r1+2*h_e,r2,sigma2));
             temp3=log(gauss(r2+2*h_e,r1,sigma2));
             temp4=log(gauss(r2-2*h_e,r1,sigma2));             
             %na kataliksei sto S1            
             cost21=tempd2+temp1+log(.5);        
             cost31=tempd3+temp1+log(.5);
             aa1=[cost21 cost31];
             cost11=min(aa1)-A;
             cost41=min(aa1)-A;
             aaa=[cost11 cost21 cost31 cost41];
             dmax(1)=max(aaa);
             [~,maxi1]=max(aaa);
             dd(1,(ii+3)/4)=maxi1; %apo pou katelikse sto S1
             % na kataliksei sto S2
             cost12=tempd1+temp2+log(.5);             
             cost42=tempd4+temp2+log(.5);
             aa1=[cost12 cost42];
             cost22=min(aa1)-A;        
             cost32=min(aa1)-A;
             aaa=[cost12 cost22 cost32 cost42];
             dmax(2)=max(aaa);
             [~,maxi2]=max(aaa);
             dd(2,(ii+3)/4)=maxi2; % apo pou katelikse sto S2
              % na kataliksei sto S3             
             cost23=tempd2+temp3+log(.5);        
             cost43=tempd4+temp3+log(.5);
             aa1=[cost23 cost43];
             cost13=min(aa1)-A;        
             cost33=min(aa1)-A;
             aaa=[cost13 cost23 cost33 cost43];
             dmax(3)=max(aaa);
             [~,maxi3]=max(aaa);
             dd(3,(ii+3)/4)=maxi3; %apo pou katelikse sto S3
              % na kataliksei sto S4
             cost14=tempd1+temp4+log(.5);                   
             cost34=tempd3+temp4+log(.5);
             aa1=[cost14 cost34];
             cost24=min(aa1)-A;        
             cost44=min(aa1)-A;
             aaa=[cost14 cost24 cost34 cost44];
             dmax(4)=max(aaa);
             [~,maxi4]=max(aaa);
             dd(4,(ii+3)/4)=maxi4; %apo pou katelikse sto s4
             %
             tempd1=dmax(1);
             tempd2=dmax(2);
             tempd3=dmax(3);
             tempd4=dmax(4);
            end
            %
        end
        %detect
        [~,maxi]=max(dmax);        
        if(maxi==1||maxi==2)
         detect(nbits)=0;
        else
         detect(nbits)=1;
        end
        for ii=nbits:-1:2
         if(dd(maxi,ii)==1||dd(maxi,ii)==2)
          detect(ii-1)=0;
         else
          detect(ii-1)=1;
         end 
         maxi=dd(maxi,ii);
        end        
       errors=0;
       for ii=1:nbits
            if(detect(ii)~=bits(ii))
             errors=errors+1;
            end
       end
       biterror=biterror+errors;
    end
    pithsymbol(i)=biterror/(sims*(nbits))
end
semilogy(snr,pithsymbol)
grid on
xlabel('snr(db)')
ylabel('Ber')
title(  'Miller-2 Viterbi')