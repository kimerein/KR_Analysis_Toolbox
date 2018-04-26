function fitout=fitDoubleGaussian(x,y,sepMin)
% x should have first column as first data trial
% different samples are different columns
% different degrees are different rows

forceEqualVariance=0;
fitOffset=1;
forceModes180degreesApart=1; % to use, others must be 0, 1
doubleFinalPoint=1;

degvals=[0; 45; 90; 135; 180; 225; 270; 315];

yForS=[];
if size(y,2)>1
    yForS=y;
    y=mean(y,2);
end

if isempty(x)
    if doubleFinalPoint==0
        x=degvals;
        newy=[y(sepMin:end); y(1:sepMin-1)];
        newyForS=[yForS(sepMin:end,:); yForS(1:sepMin-1,:)];
        y=newy;
    else
        x=[degvals; 360];
        newy=[y(sepMin:end); y(1:sepMin)];
        newyForS=[yForS(sepMin:end,:); yForS(1:sepMin,:)];
        clear y yForS
        y=newy;
    end
end

my=0;
% Remove min y
my=min(y);
y=y-my;
if ~isempty(newyForS)
    yForS=newyForS-my;
end
 
% Zero pad
% x=[-45; x];
% y=[0; y];

% dblgauss=@(a,d,b,c,x,y) a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+90))^2)/(2*c^2));
% g = fittype( @(a,b,c,d,x,y) a*x.^2+b*x+c*exp( -(y-d).^2 ), ...
%     'problem', {'c','d'}, 'independent', {'x', 'y'}, ...
%     'dependent', 'z' );
% ft=fittype( @(a,d,b,c,x,y) a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+90))^2)/(2*c^2)), ...
%     'problem', {'a','b','c','d'},'independent',{'x'},'dependent','y');
% ft=fittype('@(a,d,b,c,x,y) a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+90))^2)/(2*c^2))');
% ft=fittype('dblgauss');

% ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+90))^2)/(2*c^2))','options',fo);
if forceEqualVariance==1
    if fitOffset==0
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[0.5, 90, 45, 0.5, 270]);
        ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(e))^2)/(2*c^2))','options',fo);
    else
        if forceModes180degreesApart==0
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',[0.5, 90, 45, 0.5, 270, 0.2]);
            ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(e))^2)/(2*c^2)) + f','options',fo);
        else
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',[0.5, 90, 45, 0.5, 0.2]);
            ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+180))^2)/(2*c^2)) + e','options',fo);
        end
    end
else
    if fitOffset==0
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[0.5, 90, 45, 0.5, 270, 45]);
        ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(e))^2)/(2*f^2))','options',fo);
    else
        if forceModes180degreesApart==1
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',[0.5, 90, 45, 0.5, 45, 0.2]);
            ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(b+180))^2)/(2*e^2)) + f','options',fo);
        else
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',[0.5, 90, 45, 0.5, 270, 45, 0.2]);
            ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d*exp(-((x-(e))^2)/(2*f^2)) + g','options',fo);
        end
    end
end
[fitout,gof,output]=fit(x,y,ft);

figure(); 
hax=axes();
hl=plot(x,y+my,'Color','b');
if ~isempty(yForS)
    addErrBar(x,y+my,std(yForS,0,2)./sqrt(size(yForS,2)),'y',hax,hl);
end
hold on;
usex=linspace(x(1),x(end)+100,1000);
if forceEqualVariance==1
    if fitOffset==0
        fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.e)).^2)/(2*fitout.c^2));
    else
        if forceModes180degreesApart==0
            fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.e)).^2)/(2*fitout.c^2))+fitout.f;
        else
            fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.b+180)).^2)/(2*fitout.c^2))+fitout.e;
        end
    end
else
    if fitOffset==0
        fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.e)).^2)/(2*fitout.f^2));
    else
        if forceModes180degreesApart==1
            fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.b+180)).^2)/(2*fitout.e^2))+fitout.f;
        else
            fity=fitout.a*exp(-((usex-fitout.b).^2)/(2*fitout.c^2)) + fitout.d*exp(-((usex-(fitout.e)).^2)/(2*fitout.f^2))+fitout.g;
        end
    end
end
plot(usex,fity+my,'m');
