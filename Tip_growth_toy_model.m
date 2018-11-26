%Simple 1D (toy model) of oscillatory pollen tube growth as published in:
%Math. Model. Nat. Phenom.
%Vol. 8, No. 4, 2013, pp. 25{34
%DOI: 10.1051/mmnp/20138403
%Author: Jens Kroeger (2013)
%Free for reproduction with proper citation



%Declaration of main variables of the model.

deltat=0.0025;%time step for numerical integration
pressure=0.3; %Cytoplasmic pressure in MPa
r_0=6.0;%outter radius of pollen tube in microns
r_i=5.5;%inner radius of pollen tube in microns
caraclength=12.0;%characteristic length in microns

a_1=0.0145;%calcium absorption constant in min^{-1}
a_3=(1.5-1.5*(r_i/r_0)^2.0);
a_2=0.2;
a_6=1.5;%steepness of the sigmoid function characterizing the 
%stretch-activated calcium channels in seconds/micron
%variables for model I from 
%ROP1 is given by x
x=0.1;
%calcium is given by y
y=1;
betax=0.1;
betay=0.25;
alphax=0.25;
alphay=0.05;
xzero=0.1;
yzero=0.8;
Rzero=5;
Czero=3;
theta=8;
kappa=5;
b=0.05;


conductanceopen=150.0;%1200in picoSiemens

calcout=15000.0;%15000microMolar
calc=100.0;%microMolar
thresholdgrowthrate=63.0;
%thresholdgrowthrate=63.0/caraclength;
%micron/minute
time=256000;
thickness=zeros(1,time);
calctime=zeros(1,time);
growthrate=zeros(1,time);
viscositytime=zeros(1,time);
longtimeaxis=1:time;

xtime=zeros(1,time);
ytime=zeros(1,time);

%Beginning of loop over 4 different cases with different cell wall viscosity.

for ttt=1:4

viscositymax=0.5-0.025*ttt;
viscositynew=0.2*viscositymax;%Mpa min
viscosity=0.8*viscositymax;%Mpa min

%Beginning of loop time for solution of numerical integration. The variables
%not be recorded over the first 9000 time steps as the feedback loop stabilizes
%and the variables begin to oscillate regularly.

for i=1:time
calctime(i)=1.0;%microMolar
growthrate(i)=19.0;%microns per min
viscositytime(i)=10.0;%in MPa second
thickness(i)=1.7;%microns

end

reso=round(time/8);

thicknessresol=zeros(1,reso);
calctimeresol=zeros(1,reso);
growthrateresol=zeros(1,reso);
viscositytimeresol=zeros(1,reso);
timeaxis=zeros(1,reso);


for i=1:reso
timeaxis(i)=i*deltat*60.0;
thicknessresol(i)=1.7;
calctimereso(i)=1.0;%microMolar
growthrateresol(i)=19.0;%microns per min
viscositytimereso(i)=10.0;%in MPa second
end

%Solution of differential equations using simple Newton model. Iteration of 
%variables according to model described in Math. Model. Nat. Phenom.
%Vol. 8, No. 4, 2013, pp. 25{34
%DOI: 10.1051/mmnp/20138403

for t=90000:time
ytime(t)=y;
xtime(t)=x;
    y=y+deltat*(betay*(xtime(t-ceil(theta/60/deltat))+b)-alphay*y);
    fofx=0;
    if(x<Rzero)
        fofx=Rzero-x;
    end
    gofx=0;
        if y<Czero
       gofx=kappa*x*(Czero-y)/Rzero/Czero;
        end
        
        x=x+deltat*(betax*fofx*gofx-alphax*x);
        strainrate=growthrate(t-1)/caraclength;
        
conductance=conductanceopen/(1.0+exp(-a_6*(caraclength*strainrate-caraclength*thresholdgrowthrate/caraclength)));
R=1.5/pi*a_2*calc;%net vesicle fusion rate in microns per second
secretion=0.05/R;
thickness(t)=thickness(t-1)-deltat*secretion+0.01*(1-thickness(t-1)/viscositymax);
calc=calc+deltat*conductance*(calcout-calc)-a_1*calc;
viscosity=viscositytime(t-1)-deltat*viscositynew*R/thickness(t)+0.01*(1-viscosity/viscositymax);

%Recording of variables for plots

calctime(t)=calc;
growthrate(t)=pressure*caraclength/viscosity/thickness(t)-0.1;
viscositytime(t)=viscosity;
end
growthrate=growthrate/3.0;
for i=1:reso
thicknessresol(i)=thickness(i+time-reso);
calctimeresol(i)=calctime(i+time-reso);
growthrateresol(i)=growthrate(i+time-reso);
viscositytimeresol(i)=thickness(i+time-reso);

end
growthrateresol=growthrateresol/60.0;% in microns/s
viscositytimeresol=viscositytimeresol*60.0;% in MPa s

%Recording of arrays for figures A, B, C and D

if ttt==1
thicknessresolA=thicknessresol;
calctimeresolA=calctimeresol;
growthrateresolA=growthrateresol;
viscositytimeresolA=viscositytimeresol;
elseif ttt==2
    growthrateresolB=growthrateresol;
    viscositytimeresolB=viscositytimeresol;
elseif ttt==3
    
    growthrateresolC=growthrateresol;
    viscositytimeresolC=viscositytimeresol;
elseif ttt==4
growthrateresolD=growthrateresol;
viscositytimeresolD=viscositytimeresol;
end

%Normalization of curves 

growthrate=growthrate/60.0;
viscositytime=viscositytime*60.0;


viscositynorm=(viscositytimeresol-mean(viscositytimeresol))/max(viscositytimeresol)*5.0;
calcnorm=(calctimeresol-mean(calctimeresol))/max(calctimeresol);
thicknorm=(thicknessresol-mean(thicknessresol))/max(thicknessresol)*1.0;
growthnorm=(growthrateresol-mean(growthrateresol))/max(growthrateresol)*1.0;
end



%Output plot

figure(1);

subplot(2,2,1);
[ax, h1, h2] = plotyy(timeaxis,growthrateresolA,timeaxis,viscositytimeresolA,'plot');
title('D');
set(h1,'Linewidth',2,'color','black')
set(h2,'LineStyle','--','Linewidth',2,'color','black')
set(h2,'color','black')
axes(ax(1))
axis([0 140 0 0.5]);
axes(ax(2))
axis([0 140 0 30]);
xlabel('t (s)'); % set the X label
%ylabel('v ({\mu}m/s)'); % set the Y label
set(get(ax(1), 'Ylabel'), 'String', 'v ({\mu}m/s)')
set(get(ax(2), 'Ylabel'), 'String', '\eta (MPa s)');

subplot(2,2,2);
plot(timeaxis,growthrateresolB,'color','black','Linewidth',2);
axis([0 140 0 0.5]);
title('C');
xlabel('t (s)'); % set the X label
ylabel('v ({\mu}m/s)'); % set the Y label

subplot(2,2,3);
plot(timeaxis,growthrateresolC,'color','black','Linewidth',2);
axis([0 140 0 0.5]);
title('B');
xlabel('t (s)'); % set the X label
ylabel('v ({\mu}m/s)'); % set the Y label
 
subplot(2,2,4);
plot(timeaxis,growthrateresolD,'color','black','Linewidth',2);
title('A');
axis([0 140 0 0.5]);
xlabel('t (s)'); % set the X label
ylabel('v ({\mu}m/s)'); % set the Y label
