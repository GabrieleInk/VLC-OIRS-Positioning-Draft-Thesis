function f_scen = plotCoverage(LED, RISs, x_max, y_max, z_max, alphas, betas, Ns, ws, hs)


res=0.01;
rem=zeros(x_max/res + 1,y_max/res + 1);

f_scen=figure();
f_scen.Position= [20,20,700,700];
scatter3(LED(1),LED(2),LED(3), 'MarkerFaceColor',		[1 0 0])
hold on

for i=1:size(RISs,1)
if abs(Ns(i,2))==1
rect=[RISs(i,:)+[ws(1)/2 0 hs(1)/2];RISs(i,:)+[-ws(1)/2 0 hs(1)/2];RISs(i,:)+[-ws(1)/2 0 -hs(1)/2];RISs(i,:)+[ws(1)/2 0 -hs(1)/2]];
else
    rect=[RISs(i,:)+[ 0 ws(1)/2 hs(1)/2];RISs(i,:)+[ 0 -ws(1)/2 hs(1)/2];RISs(i,:)+[ 0 -ws(1)/2 -hs(1)/2];RISs(i,:)+[ 0 ws(1)/2 -hs(1)/2]];

end
fill3(rect(:,1),rect(:,2),rect(:,3),[0 0.4470 0.7410])
hold on
end
a=1;

for x_scan=0:res:x_max
b=1;
for y_scan=0:res:y_max
    test=[x_scan,y_scan,0];
    for i=1:size(RISs,1)
    if isCovered(LED', RISs(i,:)', test', alphas(i), betas(i),Ns(i,:),ws(i),hs(i))
        rem(b,a)=rem(b,a)+1;
    end
    end
    b=b+1;
end
a=a+1;
end
[x,y]=meshgrid(0:res:x_max,0:res:y_max);
s=pcolor(x,y,rem);
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
xlim([0 x_max])
ylim([0 y_max])
zlim([0 z_max])
daspect([1 1 1])


end