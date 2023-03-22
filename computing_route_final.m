format long

fileofell='E:\eedata\H11elip.xlsx';
fileofcircle='E:\eedata\H11cir.xlsx';
fileofgeo='E:\eedata\H11geo.xlsx';

%great cicle para
a=6377397.155;
b=6356078.963;
covt=1/180*pi;
plax=-20*covt;
longx=10*covt;
play=40*covt;
longy=120*covt;
delta_lam=acos(sin(plax).*sin(play)+cos(plax).*cos(play).*cos(longy-longx));
disp(delta_lam);
azi_circle=asin(cos(play)*sin(longy-longx)./sin(delta_lam));
disp(azi_circle);
radius=6371393;
distance_arc=radius.*delta_lam;
disp('circle arc length');
disp(distance_arc);
es=(a.^2-b.^2)./(a.^2);
eso=(b.^2)./(a.^2);
disp(eso);

%covrt parametric to geodesic coordinate
geola=zeros(3003,3);
pala=zeros(3003,1);
geolo=zeros(3003,1);
%analysis of function domain latitude is -/+(pi/2)
%longtitude is -/+(pi)
%conversion below with not be the wrong value
glax=atan(tan(plax)./sqrt(eso));
disp(glax);
disp(plax);    
glay=atan(tan(play)./sqrt(eso));
disp(glay);
disp(play);
cons=(cos(glax).*cos(longx).*cos(glay).*sin(longy)-cos(glax).*sin(longx).*cos(glay).*cos(longy));
lp=(cos(glax).*sin(longx).*eso.*sin(glay)-cos(glay).*sin(longy).*eso.*sin(glax))./cons;
mp=(-cos(glax).*cos(longx).*eso.*sin(glay)+cos(glay).*cos(longy).*eso.*sin(glax))./cons;
disp('lp');
disp(lp);
disp('mp');
disp(mp);
fun=@(x)sqrt(1./(1+eso.*((cos(x).*lp+sin(x).*mp)./eso).^2).*(1+(1+((cos(x).*lp+sin(x).*mp)./eso).^2)./(((1+eso.*((cos(x).*lp+sin(x).*mp)./eso).^2)).^2).*(lp.*sin(x)-mp.*cos(x)).^2));
disp(lp);
disp(mp);
arc_ellipse=a.*integral(fun,longx,longy);
disp('arc length');
disp(arc_ellipse);
inc_step=(longy-longx)./3000;
disp('geodesic');

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
xlswrite(fileofell,geola);

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
disp('circle geodesic');
for i = 1:3001
    xx=longx+(i-1).*inc_step;
    pos_p_arc=atan(-(cos(xx).*clp+sin(xx).*cmp));
    cgeola(i,1)=pos_p_arc;
    cgeola(i,2)=xx;
    %disp(pos_p_arc);
    %disp(xx);
    %cpala(i)=pos_p_arc;
end   
xlswrite(fileofcircle,cgeola);

%geodesic
paravertex=acos(cos(plax).*(cos(play).*sin(longy-longx)./sin(delta_lam)));
us=(a.^2./b.^2-1).*(sin(paravertex)).^2;
%find o1 using intersection
lavertex=0;
logvertex=atan(-clp./cmp);
cir_delta_lam=acos(cos(plax).*cos(abs(longx-logvertex)));
func=@(x)sqrt(1+us.*power(sin(x),2));
%geolengthofnumericalintegral
length_geo=b.*integral(func,cir_delta_lam,cir_delta_lam+delta_lam);
disp(length_geo);

%percentage difference of circle and ellipse to geodesic
disp((distance_arc-length_geo)./length_geo);
disp((arc_ellipse-length_geo)./length_geo);

%geopath
%egoe 1 for latitude 2 for longtitude
egoe=zeros(3003,2);
sinazi_depart=cos(play)*sin(longy-longx)./sin(delta_lam);
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
xlswrite(fileofgeo,egoe);

diff=zeros(3003,1);
s=0;
for i=1:3001
    diff(i)=(egoe(i,2)-geola(i,2)).^2;
    s=s+diff(i);
end
s=s/3001;
as=sqrt(s);

disp(as);