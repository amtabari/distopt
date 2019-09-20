function [objective]= TettaHandle(F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,xxx)

format bank
%% constraints

bldlupperh=    bldlupper   ;
bhdhlowerh=    bhdhlower   ;

%% Specifications

X=Xcomp;            % mole fraction
base= 1 ;                         % base component for calculation of alfa
% base : between { 1,2,c ... , c }
F ;                          % [  lbmol/h  ]
D  ;                           % [  lbmol/h  ]
B= F-D ;                          % [  lbmol/h  ]
N  = xxx(1)  ;                            % number of stages
f = xxx(2)   ;                            % feed stage
c=numel(X);                       % number of component

%%  Comment one

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fs_uipanel23==2 || fs_uipanel23==3
    feedstate={'Buuble-Point or Subcooled Liquid'};
elseif fs_uipanel23==4 || fs_uipanel23==5
    feedstate={'Dew-Pointe or Super heated Vapor'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if xxx(3)==1
    condensor={'Total-Condensor'};
else
    condensor={'Partial-Condensor'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initiation

TF   ;                         % feed temperature [  "R   ]

T(1:N)= TF ;      % temperature      [  "R   ]

%%
% Assume Vj  ( j= 1, ... ,N )

Vinitial= D*(R+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=zeros(1,N);
V(1,:)=Vinitial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:f-2
    L(j)=V(j+1)-D;
end

feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
if feed == 1
    VF=0;
else
    VF=F;
end
L(f-1)=V(f)+VF-D;
for j=f:N-1
    L(j)=V(j+1)+B;
end
L(N)=B;

%% Enthalpy functions   %%%   Change it to function @
%                             Change it to function @
%                             Change it to function @

% Vapor                           % temperature      [  "R   ]

% H1=inline(' 17000 + 30* (T-460) ');
% H2=inline(' 13000 + 20* (T-460) ');
% H3=inline(' 800   +     (T-460) ');
%
% % Liquid
% h1=inline(' 10000 + 30* (T-460) ');
% h2=inline(' 8000  + 20* (T-460) ');
% h3=inline(' 500   +     (T-460) ');

%% While Loop
for loop=1: 15
    %%   K  values
    a1=[32.71  -9.84   -25.09  ]*10^(-2);
    a2=[-9.69  67.54   102.39  ]*10^(-5);
    a3=[6.92   -37.45  -75.22  ]*10^(-8);
    a4=[-47.36 -9.07   153.84  ]*10^(-12);


    C=[ 4e3 , 8e3 , 12e3  ];
    E=[ 4.6447e3 , 4.6447e3 , 4.6447e3 ];
    for i=1:c
        for j=1:N
            k(j,i)=C(i)*exp(-E(i)/T(j));
%             k(j,i)=feval(@kk,a1(i),a2(i),a3(i),a4(i),T(j));
        end
    end
    k;
    %% Component loop , Nesting Equations
    for i=1:c

        cond=strcmpi(condensor,'Total-Condensor');
        if cond==1 % total condensor
            A(1,i)=L(1)/D;
            vPERd(2,i)=A(1,i)+1;
        else       % partial condensor
            A(1,i)=L(1)/(k(1,i)*D);
            vPERd(2,i)=A(1,i)+1;
        end

        for j=2: N-1
            A(j,i)=L(j)/(k(j,i)*V(j));
        end
        A(N,i)=B/(k(N,i)*V(N));

        for j=2: N-2
            vPERd(j+1,i)=A(j,i)*vPERd(j,i)+1;
        end

        for j=1:N-1
            lPERd(j,i)=A(j,i)*vPERd(j,i);
        end

        s(N,i)=k(N,i)*V(N)/B;
        lPERb(N-1,i)=s(N,i)+1;
        for j=f:N-1
            s(j,i)=k(j,i)*V(j)/L(j);
        end

        for j=f:N-2
            lPERb(j,i)=s(j+1,i)*lPERb(j+1,i)+1;
        end

        vPERb(f,i)=s(f,i)*lPERb(f,i);

        feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
        if feed == 1
            VF(i)=0;
            LF(i)=F*X(i);
        else
            LF(i)=0;
            VF(i)=F*X(i);
        end
        bPERd(i)= ( lPERd(f-1,i) + LF(i)/(F*X(i)) ) / ( vPERb(f,i) + VF(i)/(F*X(i)) );
        d(i)=F*X(i)/(1+bPERd(i));
        b(i)=bPERd(i)*d(i);

    end
    vPERd;
    vPERb;
    lPERd;
    lPERb;
    s;
    bPERd;
    d;
    b;
    %% Teta Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F= 100 ;
    bPERdCA=bPERd;
    vPERdCA=vPERd;
    lPERdCA=lPERd;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    teta(1)=    0 ; % initial guess
    tolerance= 1e-4 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tol=tolerance+1;
    loop=1;
    while tol > tolerance
        loop=loop+1;

        gteta(loop-1)=0;
        for i=1:c % number of component
            gteta(loop-1)=gteta(loop-1)+ F*X(i)/(1+teta(loop-1)*bPERdCA(i));
        end
        gteta(loop-1)=gteta(loop-1)-D;

        gprime(loop-1)=0;
        for i=1:c % number of component
            gprime(loop-1)=gprime(loop-1)+ bPERdCA(i)*F*X(i)/(1+teta(loop-1)*bPERdCA(i))^2;
        end
        gprime(loop-1)=-gprime(loop-1);

        teta(loop)=teta(loop-1)-gteta(loop-1)/gprime(loop-1);
        tol=abs(gteta(loop-1));

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    teta=teta(end);

    %% Corrected Compositions   &   New Temperatures
    for i=1:c
        dCO(i)=F*X(i)/(1+ teta*bPERd(i));
    end
    dCO;
    bCO=bPERdCA.*dCO;
    Xd=dCO/D;
    Xb=teta*bPERdCA.*dCO/B;

    for i=1:c
        for j=1:N
            alfa(j,i)=k(j,i)/k(j,base);
        end
    end
    alfa;

    SUMalfaXd=sum(alfa(1,:)*Xd(:));
    K(1,base)=1/SUMalfaXd;
    T(1)=E(base)/(log(C(base)/K(1,base)));
%     T(1)=(-3563.95)/( log(K(1,base))-(4.95) );

    for i=1:c
        for j=2:N-1
            SUMvPERdCAdCO(j)=sum(vPERd(j,:)*dCO(:));
            y(j,i)=vPERdCA(j,i)*dCO(i)/SUMvPERdCAdCO(j);
        end
    end


    cond=strcmpi(condensor,'Total-Condensor');
    if cond==1 % total condensor
        x(1,:)=Xd(:);
    else       % partial condensor
        y(1,:)=Xd(:);
    end

    for i=1:c
        for j=2:N-1
            SUMlPERdCAdCO(j)=sum(lPERd(j,:)*dCO(:));
            x(j,i)=lPERdCA(j,i)*dCO(i)/SUMlPERdCAdCO(j);
        end
    end
    x(N,:)=Xb;
    y(N,:)=Xb.*k(N,:);

    for j=2:N-1
        SUMyalfa(j)=sum(y(j,:)./alfa(j,:));
        K(j,base)=SUMyalfa(j);
        T(j)=E(base)/log(C(base)/K(j,base));
%         T(j)=(-3563.95)/( log(K(j,base))-(4.95) );
    end

    SUMalfaXb=sum(alfa(N,:)*Xb(:));
    K(N,base)=1/SUMalfaXb;
    T(N)=E(base)/(log(C(base)/K(N,base)));
%     T(N)=(-3563.95)/( log(K(N,base))-(4.95) );

 

    %%  Enthalpies   &   Utilities Duty
    
    c1=[-17.89  -8.48     -12.42 ];
    c2=[1.73    1.62      1.88   ]/10;
    c3=[-3.75   -1.94     -2.48  ]*10^(-4);
    e1=[44.44   61.33     71.82  ];
    e2=[501.04  588.75    658.85 ]*10^(-4);
    e3=[7.32    11.94     11.29  ]*10^(-6);
    

    for j=1:N

                [h(j,1)]=feval(@h1,T(j));
                [h(j,2)]=feval(@h2,T(j));
                [h(j,3)]=feval(@h3,T(j));
%         for i=1:3
%             [h(j,i)]=feval(@hnew,c1(i),c2(i),c3(i),T(j));
%             [H(j,i)]=feval(@HHnew,e1(i),e2(i),e3(i),T(j));
%         end

                [H(j,1)]=feval(@HH1,T(j));
                [H(j,2)]=feval(@HH2,T(j));
                [H(j,3)]=feval(@HH3,T(j));
    end
    % Vapor Feed
%     for i=1:3
%         [HF(i)]=feval(@HHnew,e1(i),e2(i),e3(i),TF);
%         [hF(i)]=feval(@hnew,c1(i),c2(i),c3(i),TF);
%     end
    [HF(1)]=feval(@HH1,TF);
    [HF(2)]=feval(@HH2,TF);
    [HF(3)]=feval(@HH3,TF);
    % Liquid Feed
    [hF(1)]=feval(@h1,TF);
    [hF(2)]=feval(@h2,TF);
    [hF(3)]=feval(@h3,TF);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cond=strcmpi(condensor,'Total-Condensor');
    if cond==1 % total condensor
        Qc=(L(1)+D)*sum((H(2,:)-h(1,:))*Xd(:));
        for j=2:f-2
            L(j)=( Qc-D*sum( (H(j+1,:)-h(1,:)).*Xd(:)' ) )/( sum( (H(j+1,:)-h(j,:)).*x(j,:) ) );
        end
        feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
        if feed == 1
            L(f-1)=( Qc-D*sum( (H(f,:)-h(1,:)).*Xd(:)' ) )/( sum( (H(f,:)-h(f-1,:)).*x(f-1,:) ) );
        else
            L(f-1)=( Qc-D*sum( (H(f,:)-h(1,:)).*Xd(:)' ) + F.* sum( (H(f,:)-HF(:)').*X(:)' ) )/( sum( (H(f,:)-h(f-1,:)).*x(f-1,:) ) );
        end
    else       % partial condensor
        Qc=L(1)*sum(((H(2,:)-h(1,:)).* Xd(:)')/k(1,:))+D*sum((H(2,:)-H(1,:))*Xd(:)) ;
        for j=2:f-2
            L(j)=( Qc-D*sum( (H(j+1,:)-H(1,:)).*Xd(:)' ) )/( sum( (H(j+1,:)-h(j,:)).*x(j,:) ) );
        end
        feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
        if feed == 1
            L(f-1)=( Qc-D*sum( (H(f,:)-H(1,:)).*Xd(:)' ) )/( sum( (H(f,:)-h(f-1,:)).*x(f-1,:) ) );
        else
            L(f-1)=( Qc-D*sum( (H(f,:)-H(1,:)).*Xd(:)' ) + F* sum( (H(f,:)-HF(:)').*X(:)' ) )/( sum( (H(f,:)-h(f-1,:)).*x(f-1,:) ) );
        end
    end
    %% Qr , Reboiler Duty
    if cond==1 % total condensor
        feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
        if feed == 1
            Qr=B*sum(h(N,:)*Xb(:)) + D*sum(h(1,:).*Xd(:)') + Qc - F*sum(hF(:).*X(:));
        else
            Qr=B*sum(h(N,:)*Xb(:)) + D*sum(h(1,:).*Xd(:)') + Qc - F*sum(HF(:).*X(:));
        end
    else       % partial condensor
        feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
        if feed == 1
            Qr=B*sum(h(N,:)*Xb(:)) + D*sum(H(1,:).*Xd(:)') + Qc - F*sum(hF(:).*X(:));
        else
            Qr=B*sum(h(N,:)*Xb(:)) + D*sum(H(1,:).*Xd(:)') + Qc - F*sum(HF(:).*X(:));
        end
    end
    Qc;
    Qr;
    for j=f: N-1
        V(j+1)=( Qr - B*sum( (h(N,:)-h(j,:)).*Xb(:)' ) )/( sum( (H(j+1,:)-h(j,:)).*y(j+1,:) ) );
    end
    %% Corrected flows
    % V(4)=L(3)+B;
    for j=1:f-2
        V(j)=L(j)+D;
    end
    feed=strcmpi(feedstate,'Buuble-Point or Subcooled Liquid');
    if feed == 1
        V(f)=L(f-1)+D;
        V(f+1)=L(f-1)+F-B;
    else
        V(f)=L(f-1)+D-F;
        V(f+1)=L(f-1)-B;
    end

    for j=f:N-1
        L(j)=V(j+1)+B;
    end

    L;
    V;
    T'-460;
    teta;
  
    
end
%%

teta;
b=b';
d=d';

V';
T=T';
T=T-460;

%%


blight=b(3);
bheavy=b(1);
dlight=d(3);
dheavy=d(1);
bldl=blight/dlight;
bhdh=bheavy/dheavy;
sl=bldl./bldlupperh;
sh=bhdhlowerh./bhdh;
objective=((N-2).^2+(sl^2+sh^2))./100;
% objective=1/(((N-2).^2+(sl^2+sh^2)));


