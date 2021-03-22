function p=parsec(n,pars)

%Parametros Aerodinamicos
rbae=pars(1);
rbai=pars(2);
abs=pars(3);
bbs=pars(4);
zbs=pars(5);
dzbs=pars(6);
xe=pars(7);
ze=pars(8);
zxxe=pars(9);
xi=pars(10);
zi=pars(11);
zxxi=pars(12);

%Parametros PARSEC
ae=zeros(1,6);
ai=zeros(1,6);
ae(1)=sqrt(2*rbae);
ai(1)=-sqrt(2*rbai);

b=([[3/4*xe^(-1/2),15/4*xe^(1/2),35/4*xe^(3/2),63/4*xe^(5/2),99/4*xe^(7/2)];...
        [3/2*xe^(1/2),5/2*xe^(3/2),7/2*xe^(5/2),9/2*xe^(7/2),11/2*xe^(9/2)];...
        [xe^(3/2),xe^(5/2),xe^(7/2),xe^(9/2),xe^(11/2)];...
        [3/2,5/2,7/2,9/2,11/2];...
        [1,1,1,1,1]])\[zxxe+1/4*ae(1)*xe^(-3/2);-1/2*ae(1)*xe^(-1/2);ze-ae(1)*xe^(1/2);tan((2*abs-bbs)/2)-1/2*ae(1);zbs+1/2*dzbs-ae(1)]; 
    
ae(2:6)=b;

b=([[3/4*xi^(-1/2),15/4*xi^(1/2),35/4*xi^(3/2),63/4*xi^(5/2),99/4*xi^(7/2)];...
        [3/2*xi^(1/2),5/2*xi^(3/2),7/2*xi^(5/2),9/2*xi^(7/2),11/2*xi^(9/2)];...
        [xi^(3/2),xi^(5/2),xi^(7/2),xi^(9/2),xi^(11/2)];...
        [3/2,5/2,7/2,9/2,11/2];...
        [1,1,1,1,1]])\[zxxi+1/4*ai(1)*xi^(-3/2);-1/2*ai(1)*xi^(-1/2);zi-ai(1)*xi^(1/2);tan((2*abs+bbs)/2)-1/2*ai(1);zbs-1/2*dzbs-ai(1)];
    
ai(2:6)=b;

%Ordenadas
t=(0:pi/2/(n-1):pi/2)';
x=1-cos(t);

%Abscisas
ze=ae(1)*x.^(1/2)+ae(2)*x.^(3/2)+ae(3)*x.^(5/2)+ae(4)*x.^(7/2)+ae(5)*x.^(9/2)+ae(6)*x.^(11/2);
zi=ai(1)*x.^(1/2)+ai(2)*x.^(3/2)+ai(3)*x.^(5/2)+ai(4)*x.^(7/2)+ai(5)*x.^(9/2)+ai(6)*x.^(11/2);

p=[flipud([x,ze]);[x(2:n),zi(2:n)]];
