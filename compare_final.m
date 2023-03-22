format long
fileofgeo='E:\eedata\startg12geo.xlsx';
a=6377397.155;
b=6356078.963;
covt=1/180*pi;
radius=6371393;
es=(a.^2-b.^2)./(a.^2);
eso=(b.^2)./(a.^2);

delta_lam=zeros(100,1);
azi_circle=zeros(100,1);
distance_arc=zeros(100,1);
length_geo=zeros(100,1);
arc_ellipse=zeros(100,1);
as=zeros(100,1);
s=zeros(100,1);
ecd=zeros(100,1);
ged=zeros(100,1);


%same starting point

dlax=0.5;
dlox=0;
dloy=90;


plax=dlax*covt;
longx=dlox*covt;

longy=dloy*covt;

eposes=0;
icre=6;

endp=zeros(100,1);

aaa=zeros(10,1);
aaa(1)=plax;
aaa(2)=longx;
aaa(3)=longy;

%different ending point
for k=1:16
    endp(k)=eposes+(k-1).*icre;
end
endp(16)=89.5;
endp(1)=0.5;

for k=1:16
    %great cicle para



play=endp(k).*covt;

delta_lam(k)=acos(sin(plax).*sin(play)+cos(plax).*cos(play).*cos(longy-longx));
azi_circle(k)=asin(cos(play)*sin(longy-longx)./sin(delta_lam(k)));
distance_arc(k)=radius.*delta_lam(k);




%covrt parametric to geodesic coordinate
geola=zeros(3003,3);
pala=zeros(3003,1);
geolo=zeros(3003,1);
%analysis of function domain latitude is -/+(pi/2)
%longtitude is -/+(pi)
%conversion below with not be the wrong value
glax=atan(tan(plax)./sqrt(eso));

glay=atan(tan(play)./sqrt(eso));

cons=(cos(glax).*cos(longx).*cos(glay).*sin(longy)-cos(glax).*sin(longx).*cos(glay).*cos(longy));
lp=(cos(glax).*sin(longx).*eso.*sin(glay)-cos(glay).*sin(longy).*eso.*sin(glax))./cons;
mp=(-cos(glax).*cos(longx).*eso.*sin(glay)+cos(glay).*cos(longy).*eso.*sin(glax))./cons;
fun=@(x)sqrt(1./(1+eso.*((cos(x).*lp+sin(x).*mp)./eso).^2).*(1+(1+((cos(x).*lp+sin(x).*mp)./eso).^2)./(((1+eso.*((cos(x).*lp+sin(x).*mp)./eso).^2)).^2).*(lp.*sin(x)-mp.*cos(x)).^2));

arc_ellipse(k)=a.*integral(fun,longx,longy);

inc_step=(longy-longx)./3000;


for i = 1:3001
    xx=longx+(i-1).*inc_step;
    pos_p_arc=atan(-(cos(xx).*lp+sin(xx).*mp)./eso);
    geola(i,1)=pos_p_arc;
    geola(i,2)=xx;
    %disp(pos_p_arc);
    %disp(xx);
    
    geola(i,3)=atan(tan(pos_p_arc).*(eso.^(0.5)));
end   
%1-geola,2-geolo,3-parala into excel table in rad


%above is the coordinate for great ellipse, below is for great circle,
%and e^2=0
cgeola=zeros(3003,3);
cpala=zeros(3000,1);
cgeolo=zeros(3000,1);
cglax=plax;
cglay=play;
clp=(cos(cglax).*sin(longx).*sin(cglay)-cos(cglay).*sin(longy).*sin(cglax))./(cos(cglax).*cos(longx).*cos(cglay).*sin(longy)-cos(cglax).*sin(longx).*cos(cglay)*cos(longy));
cmp=(-cos(cglax).*cos(longx).*sin(cglay)+cos(cglay).*cos(longy).*sin(cglax))./(cos(cglax).*cos(longx).*cos(cglay).*sin(longy)-cos(cglax).*sin(longx).*cos(cglay)*cos(longy));


inc_step=(longy-longx)./3000;
for i = 1:3001
    xx=longx+(i-1).*inc_step;
    pos_p_arc=atan(-(cos(xx).*clp+sin(xx).*cmp));
    cgeola(i,1)=pos_p_arc;
    cgeola(i,2)=xx;
    %disp(pos_p_arc);
    %disp(xx);
    %cpala(i)=pos_p_arc;
end   


%geodesic
paravertex=acos(cos(plax).*(cos(play).*sin(longy-longx)./sin(delta_lam(k))));
us=(a.^2./b.^2-1).*(sin(paravertex)).^2;
%find o1 using intersection
lavertex=0;
logvertex=atan(-clp./cmp);
cir_delta_lam=acos(cos(plax).*cos(abs(longx-logvertex)));
func=@(x)sqrt(1+us.*power(sin(x),2));
%geolengthofnumericalintegral
length_geo(k)=b.*integral(func,cir_delta_lam,cir_delta_lam+delta_lam(k));


%percentage difference of circle and ellipse to geodesic
ecd(k)=abs((distance_arc(k)-length_geo(k))./length_geo(k));
ged(k)=abs((arc_ellipse(k)-length_geo(k))./length_geo(k));

%geopath
%egoe 1 for latitude 2 for longtitude
egoe=zeros(3003,2);
sinazi_depart=cos(play)*sin(longy-longx)./sin(delta_lam(k));
cosla_highvertex=sinazi_depart*cos(plax);
sinla_highvertexsqr=1-cosla_highvertex.^2;
egoe(1,1)=cgeola(1,1);
egoe(1,2)=cgeola(1,2);
for i= 2:3001
    delta_lam1=acos(sin(plax).*sin(cgeola(i,1))+cos(plax).*cos(cgeola(i,1)).*cos(cgeola(i,2)-longx));
    clongd=cgeola(i,2)-longx;
    funcc=@(x)(sqrt(1-es*((1-(sin(x).^2).*sinla_highvertexsqr).^2))-1)./((1-(sin(x).^2).*sinla_highvertexsqr).^2);
    elongd=clongd+cosla_highvertex.*integral(funcc,cir_delta_lam,cir_delta_lam+delta_lam1);
    egoe(i,1)=cgeola(i,1);
    egoe(i,2)=elongd+longx;
end;

diff=zeros(3003,1);
s(k)=0;
for i=1:3001
    diff(i)=(egoe(i,2)-geola(i,2)).^2;
    s(k)=s(k)+diff(i);
end
s(k)=s(k)/3001;
as(k)=sqrt(s(k));


end
xlswrite(fileofgeo,as,1);
xlswrite(fileofgeo,ecd,2);
xlswrite(fileofgeo,ged,3);
xlswrite(fileofgeo,aaa,4);