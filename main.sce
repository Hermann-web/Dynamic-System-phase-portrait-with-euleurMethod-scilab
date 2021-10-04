function [liste_y,liste_t]=RK4(f,y0,N,T)
    y=y0;
    liste_y=[y0];
    t=0;
    liste_t=[0];
    h=T/N;
    for k = 1:N
    k1=h*f(t,y);
    k2=h*f(t+h/2,y+(k1/2));
    k3=h*f(t+h/2,y+(k2/2));
    k4=h*f(t+h,y+k3);
    y=y+(k1/6)+(k2/3)+(k3/3)+(k4/6);
    t=t+h;
    liste_y=[liste_y,y];
    liste_t=[liste_t,t];
end
endfunction

function yp=laser(t,y)
    yp=[(y(2)-1)*y(1);G*(A-y(2)-y(2)*y(1))];
endfunction

function [y]=pasEuler(systeme,ne,h,t,yn),
   dy = systeme(t,yn);
   for e=1:ne,
    y(e) = yn(e)+h*dy(e)
   end;
endfunction

function [mt,my]=euler(yi,T,N,systeme,ne),
   h = T/N;
   my = zeros(ne,N+1);
   mt = zeros(1,N+1);
   my(:,1)=yi;
   mt(1)=0;
   for k=1:N,
    my(:,k+1)=pasEuler(systeme,ne,h,mt(k),my(:,k));
    mt(k+1) = mt(k)+h;
   end;
endfunction

yi=[0.01;1];
ne=2;


T=4000;
N=5000;
y0=[0.01;1];
G=0.001;
A=1.1;

[mt,my] = euler(yi,T,N,laser,ne);
[sol,tps]=RK4(laser,y0,N,T);

subplot(2,2,1)
plot(tps,sol(1,:),tps,sol(2,:))
title('solutions RK4')
subplot(2,2,3)
plot(mt,my(1,:),mt,my(2,:))
title('solutions Euler')

subplot(2,2,2)
plot(sol(1,:),sol(2,:))
title('portrait de phase RK4')

subplot(2,2,4)
plot(my(1,:),my(2,:))
title('portrait de phase Euler')
