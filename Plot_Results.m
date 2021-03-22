%*************Plot Results*************************************************
%This code serves as the post-procesor. It plots the data obtained to
%generate the Pareto Fronts, airfoil shapes and the variation of
%aerodynamic efficiency vs. variation in Mach.

close all

fs = 24;

MEAN = solution(1:32,13);
STD = solution(1:32,14);

x = MEAN;
y = STD;

x = sort(x);
y = sort(y);

% x = x(find(x(:,1)<25),1);
% y = y(1:20,1);

plot(x,y,'o k')
% title('Pareto Front case-LHS','FontSize',fs, 'FontName','SansSerif')
xlabel('L/D Mean Value','FontSize',fs, 'FontName','SansSerif')
ylabel('L/D Standard Deviation','FontSize',fs, 'FontName','SansSerif')

hold on

% plot(x(32,1),y(32,1),'^','MarkerSize',20)
% plot(x(26,1),y(26,1),'d','MarkerSize',20)
% plot(x(22,1),y(22,1),'.','MarkerSize',20)

plot(x(32,1),y(32,1),'^ k','MarkerFaceColor','k','MarkerSize',8)
plot(x(26,1),y(26,1),'d k','MarkerFaceColor','k','MarkerSize',8)
plot(x(22,1),y(22,1),'. k','MarkerSize',20)


% plot(x(30,1),y(30,1),'^ k','MarkerFaceColor','k','MarkerSize',8)
% plot(x(29,1),y(29,1),'d k','MarkerFaceColor','k','MarkerSize',8)
% plot(x(24,1),y(24,1),'. k','MarkerSize',20)
% plot(x(20,1),y(20,1),'p k','MarkerFaceColor','k','MarkerSize',8)
% plot(x(16,1),y(16,1),'h k','MarkerFaceColor','k','MarkerSize',8)
% plot(x(1,1),y(1,1),'s k','MarkerFaceColor','k','MarkerSize',8)
hold on

grid on

% z(:,1) = (x(:,1)-15)/.3;
% n1=plot(x,z,'--');
% Sigma1 = '3\sigma';
% legend(n1,Sigma1);
% hold on
% axis([45 50 2.5 4])
w6 = (x(:,1)-42)/.6;
m6=plot(x,w6,'- k');

hold on

w5 = (x(:,1)-42)/.5;
m5=plot(x,w5,'-- k');

hold on

w4 = (x(:,1)-42)/.4;
m4=plot(x,w4,': k');

hold on

w3 = (x(:,1)-42)/.3;
m3=plot(x,w3,'-. k');

Sigma6 = '6\sigma';
Sigma5 = '5\sigma';
Sigma4 = '4\sigma';
Sigma3 = '3\sigma';
hold on
legend('Aerofoils','Aerofoil A', 'Aerofoil B', 'Aerofoil C', Sigma6,...
    Sigma5,Sigma4,Sigma3,'FontSize',fs, 'FontName','SansSerif')%, ...
       %'Aerofoil D', 'Aerofoil E', 'Aerofoil F', Sigma2);
% axis([8 19 0 3])
%axis([0 53 0 3.5])
%axis([30 53 1.5 4.5])
%axis([8 53 0.45 2.8])
axis([8 53 0.5 4])

a = find(solution(:,13)==x(30,1));
b = find(solution(:,13)==x(29,1));
c = find(solution(:,13)==x(24,1));
% d = find(solution(:,13)==x(20,1));
% e = find(solution(:,13)==x(16,1));
% f = find(solution(:,13)==x(1,1));

pars = solution(a,1:12);
pars2 = solution(b,1:12);
pars3 = solution(c,1:12);
% pars4 = solution(d,1:12);
% pars5 = solution(e,1:12);
% pars6 = solution(f,1:12);

coord = parsec(80, pars);
coord2 = parsec(80, pars2);
coord3 = parsec(80, pars3);
% coord4 = parsec(80, pars4);
% coord5 = parsec(80, pars5);
% coord6 = parsec(80, pars6);

figure
plot(coord(:,1),coord(:,2),'k')
aa = 'Aerofoil A'; 
% legend(aa);
hold on
plot(coord2(:,1),coord2(:,2),'-- k')
ab = 'Aerofoil B'; 
% legend(ab);
hold on
plot(coord3(:,1),coord3(:,2),'-+ k')
ac = 'Aerofoil C';
% legend(ac);
hold on
% plot(coord4(:,1),coord4(:,2),'-. k')
% ad = 'Aerofoil D';
% % legend(ac);
% hold on
% plot(coord5(:,1),coord5(:,2),'-+ k')
% ae = 'Aerofoil E';
% % legend(ac);
% hold on
% plot(coord6(:,1),coord6(:,2),'-x k')
% af = 'Aerofoil F';
% % legend(ac);
% hold on
% title('Case-LHS Aerofoils','FontSize',fs, 'FontName','SansSerif')
xlabel('X/C','FontSize',fs, 'FontName','SansSerif')
ylabel('Y/C','FontSize',fs, 'FontName','SansSerif')
axis([-0.02,1.05,-.05 0.15])
grid on
legend(aa,ab,ac,'FontSize',fs, 'FontName','SansSerif');
% ,ad,ae,af);

Mach = [0.25,0.3,0.4,0.425,0.45,0.475,0.5,0.525,0.6,0.7];

coord_eva = xfoil(coord,2.0,1e5,Mach, ...
                  'PCOP','PANE','oper iter 120',...
                               'oper/vpar n 9');
coord_eva2 = xfoil(coord2,2.0,1e5,Mach, ...
                  'PCOP','PANE','oper iter 120',...
                               'oper/vpar n 9');
coord_eva3 = xfoil(coord3,2.0,1e5,Mach, ...
                  'PCOP','PANE','oper iter 120',...
                               'oper/vpar n 9');
% coord_eva4 = xfoil(coord4,2.0,1e5,Mach, ...
%                   'PCOP','PANE','oper iter 120',...
%                                'oper/vpar n 9');
% coord_eva5 = xfoil(coord5,2.0,1e5,Mach, ...
%                   'PCOP','PANE','oper iter 120',...
%                                'oper/vpar n 9');
% coord_eva6 = xfoil(coord6,2.0,1e5,Mach, ...
%                   'PCOP','PANE','oper iter 120',...
%                                'oper/vpar n 9');                           
n1 = coord_eva.CL./coord_eva.CD;
n2 = coord_eva2.CL./coord_eva2.CD;
n3 = coord_eva3.CL./coord_eva3.CD;
% n4 = coord_eva4.CL./coord_eva4.CD;
% n5 = coord_eva5.CL./coord_eva5.CD;
% n6 = coord_eva6.CL./coord_eva6.CD;
sig = ones(size(Mach)).*42;
figure
plot(Mach,n1,'k')
aa = 'Aerofoil A';
hold on
plot(Mach,n2,'-- k')
ab = 'Aerofoil B';
hold on
plot(Mach,n3,'-x k')
ac = 'Aerofoil C';
hold on
plot(Mach,sig,'k','LineWidth',2)
hold on
% plot(Mach,n3,'-. k')
% ad = 'Aerofoil D';
% hold on
% plot(Mach,n3,'-+ k')
% ae = 'Aerofoil E';
% hold on
% plot(Mach,n3,'-x k')
% af = 'Aerofoil F';
% hold on
% title('L/D variation over M\infty case-LHS','FontSize',fs, 'FontName','SansSerif')
xlabel('M\infty','FontSize',fs, 'FontName','SansSerif')
ylabel('L/D','FontSize',fs, 'FontName','SansSerif')
legend(aa,ab,ac);
grid on



