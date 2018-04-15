function [bits, y_miller] = RN_miller(n)
y_miller=zeros(4*n,1);
bits=randi([0,1],n,1);
a=rand();
if(bits(1)==0)
 if(a>=.5)  
  y_miller(1)=1;
  y_miller(2)=-1;
  y_miller(3)=1;
  y_miller(4)=-1;
 elseif(a<.5)
  y_miller(1)=-1;
  y_miller(2)=1;
  y_miller(3)=-1;
  y_miller(4)=1;
 end 
elseif(bits(1)==1)
 if(a>=.5)   
  y_miller(1)=1;
  y_miller(2)=-1;
  y_miller(3)=-1;
  y_miller(4)=1;
 elseif(a<.5)
  y_miller(1)=-1;
  y_miller(2)=1;
  y_miller(3)=1;
  y_miller(4)=-1;
 end
end
for i=2:n
 if(bits(i-1)==1)
  if(bits(i)==0&&y_miller(4*(i-1)-3)==1)
   y_miller(4*i-3)=-1;
   y_miller(4*i-2)=1;
   y_miller(4*i-1)=-1;
   y_miller(4*i)=1;
  elseif(bits(i)==0&&y_miller(4*(i-1)-3)==-1)   
   y_miller(4*i-3)=1;
   y_miller(4*i-2)=-1;
   y_miller(4*i-1)=1;
   y_miller(4*i)=-1;
  elseif(bits(i)==1&&y_miller(4*(i-1)-3)==1)   
   y_miller(4*i-3)=-1;
   y_miller(4*i-2)=1;
   y_miller(4*i-1)=1;
   y_miller(4*i)=-1;
  elseif(bits(i)==1&&y_miller(4*(i-1)-3)==-1)   
   y_miller(4*i-3)=1;
   y_miller(4*i-2)=-1;
   y_miller(4*i-1)=-1;
   y_miller(4*i)=1; 
  end   
 elseif(bits(i-1)==0)
  if(bits(i)==0&&y_miller(4*(i-1)-3)==1)
   y_miller(4*i-3)=-1;
   y_miller(4*i-2)=1;
   y_miller(4*i-1)=-1;
   y_miller(4*i)=1;
  elseif(bits(i)==0&&y_miller(4*(i-1)-3)==-1)   
   y_miller(4*i-3)=1;
   y_miller(4*i-2)=-1;
   y_miller(4*i-1)=1;
   y_miller(4*i)=-1;
  elseif(bits(i)==1&&y_miller(4*(i-1)-3)==1)   
   y_miller(4*i-3)=1;
   y_miller(4*i-2)=-1;
   y_miller(4*i-1)=-1;
   y_miller(4*i)=1;
  elseif(bits(i)==1&&y_miller(4*(i-1)-3)==-1)   
   y_miller(4*i-3)=-1;
   y_miller(4*i-2)=1;
   y_miller(4*i-1)=1;
   y_miller(4*i)=-1; 
  end
  %
 end
 %
end
end