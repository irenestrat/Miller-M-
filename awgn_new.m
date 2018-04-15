clc
clear all
close all
L=100;
snr=10;
sims=20000;
nbits=16;
pithsymbol=zeros(length(snr),1);
for i=1:length(snr)
    biterror=0;    
    for j=1:sims        
        [bits,y_miller] = RN_miller(nbits);    %dimiourgia 20 bits kai ta kolame meta to preamble
        x = [y_miller];
        x_upsample = upsample(x,L/2);
        %metro=1;
        sigma = sqrt(L/(10^(snr(i)/10)));
        y_conv = conv(x_upsample(1:end),ones(L/2,1));
        n = sigma*randn(size(y_conv)) + 1i*sigma*randn(size(y_conv));
        y_r = y_conv + n;
        y_r_filt=conv(y_r,ones(L/2,1));
        y_r_sampled  = y_r_filt(L/2:L/2:length(y_conv));
        y_info=y_r_sampled(1:end);
        %tora to t*
        %tora tha kanoume detect
        detect=zeros(nbits,1);
        for ii=1:4:4*nbits-3
            y0=y_info(ii);
            y1=y_info(ii+1);
            y2=y_info(ii+2);
            y3=y_info(ii+3);
            %
             %2 ta prosarmosmena filtra
            r1=(y0-y1+y2-y3)/(4*L);
            r2=(y0-y1-y2+y3)/(4*L);
            %
            a1=abs(r1-2)^2<=abs(r1+2)^2;
            a2=abs(r1-2)^2+abs(r2)^2<=abs(r1)^2+abs(r2-2)^2;
            a3=abs(r1-2)^2+abs(r2)^2<=abs(r1)^2+abs(r2+2)^2;
            %
            b1=abs(r1+2)^2<=abs(r1-2)^2;
            b2=abs(r1+2)^2+abs(r2)^2<=abs(r1)^2+abs(r2-2)^2;
            b3=abs(r1+2)^2+abs(r2)^2<=abs(r1)^2+abs(r2+2)^2;
            %
            c1=abs(r2-2)^2+abs(r1)^2<=abs(r2)^2+abs(r1-2)^2;
            c2=abs(r2-2)^2+abs(r1)^2<=abs(r2)^2+abs(r1+2)^2;
            c3=abs(r2-2)^2<=abs(r2+2)^2;
            %
            d1=abs(r2+2)^2+abs(r1)^2<=abs(r2)^2+abs(r1-2)^2;
            d2=abs(r2+2)^2+abs(r1)^2<=abs(r2)^2+abs(r1+2)^2;
            d3=abs(r2+2)^2<=abs(r2-2)^2;
            %
            if((d1&&d2&&d3)||(c1&&c2&&c3))
             detect((ii+3)/4)=1;
            elseif((b1&&b2&&b3)||(a1&&a2&&a3))
             detect((ii+3)/4)=0; 
            end 
        end
       errors=0;
       for ii=1:nbits
            if(detect(ii)~=bits(ii))
             errors=errors+1;
            end
       end
       biterror=biterror+errors/nbits;
    end
    pithsymbol(i)=biterror/sims
end
semilogy(snr,pithsymbol)
grid on
xlabel('snr(db)')
ylabel('Ber')
title(  'AWGN')