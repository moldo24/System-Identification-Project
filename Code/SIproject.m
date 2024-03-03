t = Moldovan(:,1);
u = Moldovan(:,2);
y = Moldovan(:,3);
x = Moldovan(:,4);

%Plot intrare/iesire y
figure
plot(t, u, t, y)
xlabel('Time')
ylabel('Amplitude')
legend('Intrare', 'Iesire fara zero' )
title('Input and Output Plot')

%Plot intrare/iesire x
figure
plot(t, u, t , x)
xlabel('Time')
ylabel('Amplitude')
legend('Intrare', 'Iesire cu un zero')
title('Input and Output Plot')



%%
ymin=180; %minim iesire 
ymax=187; %maxim iesire  
umin=178; %minim intrare 
umax=185; %maxim intrare  

Mr=(y(ymin)-y(ymax))/(u(umin)-u(umax)) %modulul de rezonanta
Tr=(t(ymax)-t(ymin))*2 %perioada de rezonanta
zeta=(sqrt((Mr-sqrt(Mr^2-1)))/2*Mr)%factorul de amortizare 
wr = (2*pi/Tr) %pulsatia de rezonanta
wn = wr/sqrt(1-2*zeta^2) %pulsatia naturala
K=mean(y)/mean(u); %=1.0081
%%
Hj=tf([K*wn^2],[1 2*zeta*wn wn^2])
%%

num=K*wn^2;
den=[1 2*zeta*wn wn^2];
Hdes=tf(num,den) 

figure
bode(Hdes)


%Conditii initiale nenule=>este necesar modelul de tip spatiul starilor

A=[0 1;-wn^2 -2*zeta*wn];
B=[0;K*wn^2];
C=[1 0];
D=0;
%convertire in spatiul starilor
sys=ss(A,B,C,D);
%generare semnal simulat
ysim=lsim(sys, u , t,[y(1),0]);

%SS:ne alegem var de strare=x1-ies x2-deriv=>mat 2/2(a)

figure
plot (t,ysim,'b',t,y,'r');
xlabel('Time')
ylabel('Amplitude')
legend('Simulated Output', 'Actual Output')
title('System Response')
hold on


J=1/sqrt(length(t))*norm(y-ysim)
empn=norm(y-ysim)/norm(y-mean(y))*100%eroare medie patratica normalizata

%%

Te=t(2)-t(1);%perioada esantionare
date=iddata(y, u, Te); %dt-ti
%verificare ordin
n4sid(date,1:10)
%% Metoda arx-auto

Marx=arx(date, [2,2,0]);%na,nb,timp mort
%2.2.0: 2-ord a,2-ord b,0-ord timp mort


figure,resid(date,Marx),shg,title('resid Marx')
figure,compare(date,Marx),shg


%% Metoda vi-inter
Mvi=iv4(date, [2 2 1]);

figure,resid(date,Mvi),shg,title('resid Mvi')
figure,compare(date,Mvi),shg

%% Metoda oe-inter ,er la ies
Moe=oe(date, [2 2 1])

figure,resid(date,Moe),shg,title('resid Moe')
figure,compare(date,Moe),shg

Hz_oe1 = tf(Moe.B,Moe.F,Te,'variable','z^-1')%functia de transfer in discret
Hs_oe1 = d2c(Hz_oe1,'zoh')%functia de transfer in continu
 



%% Metoda armax-auto

Marmax=armax(date, [2 2 2 1])

figure,resid(date, Marmax),shg,title('resid Marmax')
figure,compare(date,Marmax),shg

Hz_arx1 = tf(Marmax.B,Marmax.A,Te,'variable','z^-1') %functia de transfer in discret
Hs_arx1 = d2c(Hz_arx1,'zoh') %functia de transfer in continu



