function timi = gauss(x1,x2,s)
%timi=(1/(pi*s^2)^2)*exp((-(abs(x1))^2-(abs(x2))^2)/(s^2));
timi=exp((-(abs(x1))^2-(abs(x2))^2)/(s^2));
end
