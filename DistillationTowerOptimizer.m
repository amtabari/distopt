function varargout = DistillationTowerOptimizer(varargin)
%%% Initialization

clc

%%%%%%%%%%%%%%%%%%%%%%%%%      dp_uipanel12 Gui
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DistillationTowerOptimizer_OpeningFcn, ...
    'gui_OutputFcn',  @DistillationTowerOptimizer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before DistillationTowerOptimizer is made visible.
function DistillationTowerOptimizer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = DistillationTowerOptimizer_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    tc=get(handles.tag,'value');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   set(handles.tag, 'String', Tn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   cp=str2double(get(handles.cp,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


set(handles.figure1, 'Visible', 'on');
set(handles.figure1, 'Position', [95,25,65,15.5]);
set(handles.uipanel1, 'Visible', 'on');
set(handles.uipanel1, 'Position', [2.4,0.692,60.6,14.231]);
pause(.5)
for i=1:13
    pause(0.07)
    set(handles.text17, 'FontSize', i+1);
end
pause(.5)
set(handles.text18, 'Visible', 'on');
for i=1:8
    pause(0.1)
    set(handles.text18, 'FontSize', i+1);
end
pause(0.5)
for i=5:-1:1
    set(handles.watch, 'String', i);
    pause(.4)
    set(handles.watch, 'FontSize', 10);
    pause(.1)
    set(handles.watch, 'FontSize', 12);
    pause(.1)
    set(handles.watch, 'FontSize', 14);
    pause(.1)
    set(handles.watch, 'FontSize', 16);
    pause(.1)
    set(handles.watch, 'FontSize', 14);
    pause(.1)
    set(handles.watch, 'FontSize', 12);
    pause(.1)
    set(handles.watch, 'FontSize', 10);
    pause(.1)
    set(handles.watch, 'FontSize', 8);
%     pause(.1)
end

set(handles.uipanel1, 'Visible', 'off');
set(handles.uipanel5, 'Visible', 'on');
set(handles.figure1, 'Name', 'Process Optimizer');
set(handles.uipanel5, 'Position', [2.4,0.692,60.6,14.231]);

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%         next_uipanel5 Button ,,, uipanel5           %%%%%%%%

function next_uipanel5_Callback(hObject, eventdata, handles)

dt=get(handles.dt,'Value');
fo=get(handles.fo,'Value');
popupmenu1=get(handles.popupmenu1,'Value');

if dt==1
    if popupmenu1==1 % Distillation Tower's Optimization
        set(handles.uipanel5, 'Visible', 'off');
        set(handles.figure1, 'Name', 'Distillation Tower''s Optimization');
        set(handles.figure1, 'Position', [65,3,117.3,55.8]);
        set(handles.uipanel7, 'Visible', 'on');
        set(handles.uipanel7, 'Position', [1.2,.385,116,55.5]);
        set(handles.uipanel3, 'Visible', 'on');
        set(handles.uipanel4, 'Visible', 'on');
        set(handles.uipanel23, 'Visible', 'on');
        %-----------------------------------------------------------
        set(handles.uipanel16, 'Visible', 'off');
        set(handles.uipanel13, 'Visible', 'off');
        set(handles.uipanel17, 'Visible', 'off');
        set(handles.back_uipanel7, 'Visible', 'off');
        set(handles.uipanel26, 'Visible', 'off');
        set(handles.uipanel2, 'Position', [1.4,44.462,28,9.846]);
        set(handles.uipanel8, 'Visible', 'on');
        set(handles.uipanel8, 'Position', [1.6,28.462,28,5.308]);
        

        set(handles.text29, 'String', 'Minimum Trays');
        set(handles.text46, 'String', 'Minimum Trays');
        set(handles.text56, 'String', 'Minimum Trays');
        set(handles.text30, 'String', 'Maximum Trays');
        set(handles.text47, 'String', 'Maximum Trays');
        set(handles.text57, 'String', 'Maximum Trays');

    elseif popupmenu1==2 % Distillation Tower's Simulation
        set(handles.uipanel5, 'Visible', 'off');
        set(handles.figure1, 'Name', 'Distillation Tower''s Simulation');
        set(handles.uipanel6, 'Visible', 'on');
        set(handles.uipanel6, 'Position', [2.4,0.692,60.6,14.231]);

    end
end

if fo==1 % Function Optimization
    set(handles.uipanel5, 'Visible', 'off');
    set(handles.figure1, 'Name', 'Function''s Optimization');
    set(handles.figure1, 'Position', [50,15,157.80,28]);
    set(handles.uipanel7, 'Visible', 'on');
    set(handles.uipanel7, 'Position', [1.2,.385,156,27.7]);
    set(handles.uipanel3, 'Visible', 'off');
    set(handles.uipanel4, 'Visible', 'off');
    set(handles.uipanel23, 'Visible', 'off');
    %-------------------------------------------------------------------
    set(handles.uipanel2, 'Visible', 'on');
    set(handles.uipanel2, 'Position', [1.2,16.8,28,9.846]);
    set(handles.back_uipanel7, 'Visible', 'on');
    set(handles.back_uipanel7, 'Position', [2.3,1.2,28,2]);
    set(handles.uipanel8, 'Visible', 'off');
    set(handles.uipanel26, 'Visible', 'on');
    set(handles.uipanel26, 'Position', [2.4,3.7,28,5.308]);
    
    set(handles.uipanel16, 'Visible', 'on');
    set(handles.uipanel16, 'Position', [75,17,80.5,10]);
    set(handles.uipanel13, 'Visible', 'on');
    set(handles.uipanel13, 'Position', [75,1.24,49.8,15.615]);
    set(handles.uipanel17, 'Visible', 'on');
    set(handles.uipanel17, 'Position', [127.6,1.24,27.8,15.615]);
    %-------------------------------------------------------------------
    set(handles.text29, 'String', 'Minimum Value');
    set(handles.text46, 'String', 'Minimum Value');
    set(handles.text56, 'String', 'Minimum Value');
    set(handles.text30, 'String', 'Maximum Value');
    set(handles.text47, 'String', 'Maximum Value');
    set(handles.text57, 'String', 'Maximum Value');

end

%%%%%%%%%%%%%%%%%%%%%%         next_uipanel6 Button ,,, uipanel6           %%%%%%%%

function next_uipanel6_Callback(hObject, eventdata, handles)

set(handles.feedpanel, 'Visible', 'off');
set(handles.molefractionpanel, 'Visible', 'off');


mtm=get(handles.mtm, 'Value');
tm=get(handles.tm, 'Value');

if mtm==1

set(handles.uipanel6, 'Visible', 'off');
set(handles.figure1, 'Name', 'Mc-Cabe Thiele Method');
set(handles.figure1, 'Position', [55,15,146.28,28.3]);
set(handles.uipanel12, 'Visible', 'on');
set(handles.uipanel12, 'Position', [.8,.385,145,28.077]);
set(handles.fs, 'Value', 1);

elseif tm==1
    
    set(handles.uipanel6, 'Visible', 'off');
    set(handles.figure1, 'Name', 'Theta Method');
    set(handles.figure1, 'Position', [35,15,191.28,28.3]);
    set(handles.uipanel19, 'Visible', 'on');
    set(handles.uipanel19, 'Position', [.8,.385,190,28.077]);
    set(handles.uipanel13, 'Visible', 'off');

end



%%%%%%%%%%%%%%%%%%%%%%         back_uipanel23 Button ,            %%%%%%%%

function back_uipanel23_Callback(hObject, eventdata, handles)

set(handles.figure1, 'Position', [95,25,65,15.5]);
set(handles.uipanel7, 'Visible', 'off');
set(handles.uipanel5, 'Visible', 'on');
set(handles.uipanel5, 'Position', [2.4,0.692,60.6,14.231]);
set(handles.uipanel9, 'Visible', 'off');
set(handles.uipanel10, 'Visible', 'off');
set(handles.uipanel11, 'Visible', 'off');
set(handles.figure1, 'Name','Process Optimizer');
set(handles.uipanel13, 'Visible', 'off');
set(handles.uipanel16, 'Visible', 'off');
set(handles.uipanel17, 'Visible', 'off');
% set(handles.back_uipanel23, 'Visible', 'off');



%%%%%%%%%%%%%%%%%%%%%%         back_uipanel2  Button ,            %%%%%%%%

function back_uipanel12_Callback(hObject, eventdata, handles)

set(handles.figure1, 'Position', [95,25,65,15.5]);
set(handles.uipanel12, 'Visible', 'off');
set(handles.uipanel6, 'Visible', 'on');
set(handles.uipanel6, 'Position', [2.4,0.692,60.6,14.231]);
set(handles.figure1, 'Name','Distillation Tower''s Simulation');

%%%%%%%%%%%%%%%%%%%%%%         back_uipanel9  Button ,            %%%%%%%%

function back_uipanel19_Callback(hObject, eventdata, handles)

set(handles.figure1, 'Position', [95,25,65,15.5]);
set(handles.uipanel19, 'Visible', 'off');
set(handles.uipanel6, 'Visible', 'on');
set(handles.uipanel6, 'Position', [2.4,0.692,60.6,14.231]);
set(handles.figure1, 'Name','Distillation Tower''s Simulation');


%%%%%%%%%%%%%%%%%%%%%%         back_uipanel6  Button ,            %%%%%%%%

function back_uipanel6_Callback(hObject, eventdata, handles)

set(handles.uipanel6, 'Visible', 'off');
set(handles.uipanel5, 'Visible', 'on');
set(handles.figure1, 'Name','Process Optimizer');


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


% function reset_uipanel12_Callback(hObject, eventdata, handles)


function run_uipanel7_Callback(hObject, eventdata, handles)

ra=get(handles.ra,'Value');
ga=get(handles.ga,'Value');
pso=get(handles.pso,'Value');

F=str2double(get(handles.fr_uipanel23,'String'));
D=str2double(get(handles.d_uipanel23,'String'));
R=str2double(get(handles.r_uipanel23,'String'));
TF=str2double(get(handles.ftemp_uipanel23,'String'));
bldlupper=str2double(get(handles.bldlupper,'String'));
bhdhlower=str2double(get(handles.bhdhlower,'String'));


mfoc(1)=str2double(get(handles.mfoc11, 'String'));
mfoc(2)=str2double(get(handles.mfoc22, 'String'));
mfoc(3)=str2double(get(handles.mfoc33, 'String'));
mfoc(4)=str2double(get(handles.mfoc44, 'String'));
mfoc(5)=str2double(get(handles.mfoc55, 'String'));
mfoc(6)=str2double(get(handles.mfoc66, 'String'));
mfoc(7)=str2double(get(handles.mfoc77, 'String'));
mfoc(8)=str2double(get(handles.mfoc88, 'String'));
mfoc(9)=str2double(get(handles.mfoc99, 'String'));
mfoc(10)=str2double(get(handles.mfoc1010, 'String'));

noc_uipanel23=get(handles.noc_uipanel23, 'Value');
for i=1:noc_uipanel23
    Xcomp(i)=mfoc(i);
end
fs_uipanel23=get(handles.fs_uipanel23,'Value');



if ra==1

    % parameters
    nov_uipanel9=str2double(get(handles.nov_uipanel9,'String'));
    pp_uipanel9=str2double(get(handles.pp_uipanel9,'String'));
    nop_uipanel9=str2double(get(handles.nop_uipanel9,'String'));
    noi_uipanel9=str2double(get(handles.noi_uipanel9,'String'));
    norp_uipanel9=str2double(get(handles.norp_uipanel9,'String'));
    nord_uipanel9=str2double(get(handles.nord_uipanel9,'String'));
    mint_uipanel9=str2double(get(handles.mint_uipanel9,'String'));
    maxt_uipanel9=str2double(get(handles.maxt_uipanel9,'String'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.pw, 'Visible', 'on');
    


    pause(.001)
    %     clc
    format short
    tic
    %-------------------------------------------------------------------------
    condition={'real','real','real','real','real'}; % type of variable(s)
    %-------------------------------------------------------------------------
    varb=nov_uipanel9;%   3     ;              % number of variable(s)
    nc=     1     ;              % number of combination between   1  to (varb-2)
    %-------------------------------------------------------------------------
    spp=pp_uipanel9;%    50   ;               % size of primary population (( More than 10 ))
    %-------------------------------------------------------------------------
    np(1:varb)=nop_uipanel9;%   11   ;         % number of parts for real numbers
    npfirst=np;
    secondp =noi_uipanel9;%     1   ;          % number of intervals in " axnew "
    %-------------------------------------------------------------------------
    randomno=  1    ;            % number of steps in Randomizers
    %-------------------------------------------------------------------------
    NOP=norp_uipanel9;%    np(1) ;              % number of parts in Rain Algorithm(Rain Parts)
    N=nord_uipanel9;%      9     ;              % number Of Raindrop
    rp=     0     ;              % number of randoms
    %-------------------------------------------------------------------------
    lowerbound=[    mint_uipanel9 , mint_uipanel9  , 1   ];    % lower bound of variable(s)
    upperbound=[    maxt_uipanel9 , maxt_uipanel9  , 2   ];    % upper bound of variable(s)
    %-------------------------------------------------------------------------
    e=1 ;                        % function evaluation
    %-------------------------------------------------------------------------
    for i=1:varb
        realcond=strcmpi(condition(i),'real');
        if realcond == 1
            answer=linspace(lowerbound(i),upperbound(i),np(i));
            ax(i,1:numel(answer))=answer;
        end
        integercond=strcmpi(condition(i),'integer');
        if integercond == 1
            answer=ceil(lowerbound(i)):1:floor(upperbound(i));
            ax(i,1:numel(answer))=answer;
            np(i)=numel(answer);
        end
    end
    ax;
    axfirst=ax;


    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %                  MatrixAlgorithmVersion3 Func : MAVersion3
    [bestfxnew,bestmxnew,bestxnew,ax,np,e]=feval(@MAVersion3FuncOpt,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,...
    ax,spp,varb,secondp,np,e,nc,randomno,NOP,N,rp,lowerbound,upperbound);
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    bestfxnew';
    bestmxnew;
    bestxnew(:,end)';

    bestfxfirst=bestfxnew(end);
    bestmxfirst=bestmxnew(:,end);
    bestxfirst=bestxnew(:,end)'
    e;


    seconds=toc;

    set(handles.pw, 'Visible', 'off');

    % Results
    set(handles.not, 'String',bestxfirst(1) );
    set(handles.ft, 'String', bestxfirst(2)  );
    if bestxfirst(3)==1
        set(handles.ct, 'String',  'Total');
    else
        set(handles.ct, 'String',  'Partial');
    end
    set(handles.tc, 'String',  bestfxfirst);
    set(handles.fe, 'String', e  );
    set(handles.t, 'String', seconds  );


elseif ga==1

    % parameters
    nov_uipanel10=str2double(get(handles.nov_uipanel10,'String'));
    nop_uipanel10=str2double(get(handles.nop_uipanel10,'String'));
    nog_uipanel10=str2double(get(handles.nog_uipanel10,'String'));
    nomcg_uipanel10=str2double(get(handles.nomcg_uipanel10,'String'));
    nomcr_uipanel10=str2double(get(handles.nomcr_uipanel10,'String'));
    noec_uipanel10=str2double(get(handles.noec_uipanel10,'String'));
    mint_uipanel10=str2double(get(handles.mint_uipanel10,'String'));
    maxt_uipanel10=str2double(get(handles.maxt_uipanel10,'String'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.pw, 'Visible', 'on');
    pause(.001)

    tic
    % clc
    % figure(1)
    % clf
    %     clear all
    format short

    %------------------------        parameters        ------------------------
    % befor using this function you must specified your function in fun00.m
    % file in current directory and then set the parameters
    var=nov_uipanel10;%3;            % Number of variables (this item must be equal to the
    %   number of variables that is used in the function in
    %   fun00.m file)
    n=nop_uipanel10;%100;            % Number of population

    m0=nog_uipanel10;%20;            % Number of generations that max Value remains constant
    %   (use for termination criteria)
    nmutationG=nomcg_uipanel10;%20;                  %number of mutation children(Gaussian)
    nmutationR=nomcr_uipanel10;%20;                  %number of mutation children(random)
    nelit=noec_uipanel10;%2;                        %number of elitism children
    valuemin=ones(1,var)*   mint_uipanel10 ;      % min possible value of variables
    valuemax=ones(1,var)*  maxt_uipanel10 ;     % max possible value of variables

    %-------------------------------------------------------------------------
    nmutation=nmutationG+nmutationR;
    sigma=(valuemax-valuemin)/10;    %Parameter that related to Gaussian
    %   function and used in mutation step
    max1=zeros(nelit,var);
    parent=zeros(n,var);
    cu=[valuemin(1) valuemax(1) valuemin(2) valuemax(2)];
    for l=1:var
        p(:,l)=valuemin(l)+rand(n,1).*(valuemax(l)-valuemin(l));
    end
    initial=p;
    m=m0;
    maxvalue=ones(m,1)*-1e10;
    maxvalue(m)=-1e5;
    g=0;
    meanvalue(m)=0;
    %-------------   ****    termination criteria   ****-------------
    while abs(maxvalue(m)-maxvalue(m-(m0-1)))>0.001*maxvalue(m) &...
            (abs(maxvalue(m))>1e-10 & abs(maxvalue(m-(m0-1)))>1e-10)...
            & m<10000 & abs(maxvalue(m)-meanvalue(m))>1e-5 | m<20
        sigma=sigma./(1.05);% reducing the sigma value
        %  ------  **** % reducing the number of mutation()random   **** ----
        g=g+1;
        if g>10 & nmutationR>0
            g=0;
            nmutationR=nmutationR-1;
            nmutation=nmutationG+nmutationR;
        end

        %-------------   ****    function evaluation   ****-------------
        [p]=feval(@xModifiedGA,p);
        for i=1:n
            %         y(i)=fun00(password(i,:));
            y(i)=100/(TettaHandle(F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,p(i,:)));
            %               y(i)=(TettaHandle(password(i,:)));
        end
        s=sort(y);
        maxvalue(1:nelit)=s(n:-1:n-nelit+1);
        if nelit==0
            maxvalue(1)=s(n);
            for i=1:n
                if y(i)==maxvalue(1)
                    max1(1,:)=p(i,:);
                end
            end
        end
        for k=1:nelit
            for i=1:n
                if y(i)==maxvalue(k)
                    max1(k,:)=p(i,:);
                end
            end
        end

        y=y-min(y)*1.02;
        sumd=y./sum(y);
        meanvalue=y./(sum(y)/n);


        %-------------   ****   Selection: Rolette wheel   ****-------------
        for l=1:n
            sel=rand;
            sumds=0;
            j=1;
            while sumds<sel
                sumds=sumds+sumd(j);
                j=j+1;
            end
            parent(l,:)=p(j-1,:);
        end
        p=zeros(n,var);

        %-------------   ****    regeneration   ****-------------
        for l=1:var


            %-------------   ****    cross-over   ****-------------
            for j=1:ceil((n-nmutation-nelit)/2)
                t=rand*1.5-0.25;
                p(2*j-1,l)=t*parent(2*j-1,l)+(1-t)*parent(2*j,l);
                p(2*j,l)=t*parent(2*j,l)+(1-t)*parent(2*j-1,l);
            end


            %-------------   ****    elitism   ****-------------
            for k=1:nelit
                p((n-nmutation-k+1),l)=max1(k,l);
            end


            %-------------   ****    mutation   ****-------------
            for i=n-nmutation+1:n-nmutationR
                phi=1-2*rand;
                z=erfinv(phi)*(2^0.5);
                p(i,l)=z*sigma(l)+parent(i,l);

            end
            for i=n-nmutationR+1:n
                p(i,1:var)=valuemin(1:var)+rand(1,var).*(valuemax(1:var)-...
                    valuemin(1:var));
            end
            for i=1:n
                for l=1:var
                    if p(i,l)<valuemin(l)
                        p(i,l)=valuemin(l);
                    elseif p(i,l)>valuemax(l)
                        p(i,l)=valuemax(l);
                    end
                end
            end
        end
        p;
        m=m+1;
        max1;
        maxvalue(m)=maxvalue(1);
        maxvalue00(m-m0)=maxvalue(1);(1/maxvalue(1)).*100;
        mean00(m-m0)=sum(s)/n;
        meanvalue(m)=mean00(m-m0);

    end


    set(handles.pw, 'Visible', 'off');

%     disp('**************************************')
    num_of_fun_evaluation=n*m;
    max_point_GA=max1(1,:);
    maxvalue_GA=100/maxvalue00(m-m0);
    %     maxvalue_GA=maxvalue00(m-m0);
    seconds=toc;

    % Results
    set(handles.not, 'String', max_point_GA(1) );
    set(handles.ft, 'String',max_point_GA(2)  );
    if max_point_GA(3)==1
        set(handles.ct, 'String',  'Total');
    else
        set(handles.ct, 'String',  'Partial');
    end
    set(handles.tc, 'String',maxvalue_GA);
    set(handles.fe, 'String',num_of_fun_evaluation );
    set(handles.t, 'String', seconds  );





elseif pso==1

    set(handles.pw, 'Visible', 'on');
    pause(.001)

    % parameters
    nov_uipanel11=str2double(get(handles.nov_uipanel11,'String'));
    pa_uipanel11=(get(handles.pa_uipanel11,'String'));
    pb_uipanel11=(get(handles.pb_uipanel11,'String'));
    pc_uipanel11=str2double(get(handles.pc_uipanel11,'String'));
    pd_uipanel11=str2double(get(handles.pd_uipanel11,'String'));
    ivov_uipanel11=str2double(get(handles.ivov_uipanel11,'String'));
    ivox_uipanel11=str2double(get(handles.ivox_uipanel11,'String'));
    noi_uipanel11=str2double(get(handles.noi_uipanel11,'String'));
    mint_uipanel11=str2double(get(handles.mint_uipanel11,'String'));
    maxt_uipanel11=str2double(get(handles.maxt_uipanel11,'String'));


    parametra=strcmpi(pa_uipanel11,'set');
    parametrb=strcmpi(pb_uipanel11,'set');
    if parametra==1 && parametrb==1
        a=[.9 .7 .9 .1 .1 -.7 ];
        b=[.1 .3  3 .1 2.1 .5 ];
    else
        pa_uipanel11=str2double(get(handles.pa_uipanel11,'String'));
        pb_uipanel11=str2double(get(handles.pb_uipanel11,'String'));
        a=pa_uipanel11;%[.9 .7 .9 .1 .1 -.7 ];
        b=pb_uipanel11;%[.1 .3  3 .1 2.1 .5 ];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tic
    e=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varb=nov_uipanel11;%3;

    c=pc_uipanel11;%;
    d=pd_uipanel11;%1;

    lowerbound=mint_uipanel11;%   4   ;    % lower bound of variable(s)
    upperbound= maxt_uipanel11;%  100  ;    % upper bound of variable(s)

    for i=1:numel(a)

        v(1:varb,1)=ivov_uipanel11;%-0.1;
        x(1:varb,1)=ivox_uipanel11;%10;
        x(1:varb,1)=feval(@xModifiedPSO,x(1:varb,1),lowerbound,upperbound);
        p(1)=feval(@TettaHandle,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,x(1:varb,1));

        v(1:varb,2)=a(i)*v(1:varb,1)+b(i)*(p(1)-x(1:varb,1));
        x(1:varb,2)=c*x(1:varb,1)+d*v(1:varb,1);

        k=1;
        for j=1:noi_uipanel11;%500
            k=k+1;
            v(1:varb,k+1)=a(i)*v(1:varb,k)+b(i)*rand(varb,1).*(p(k-1)-x(1:varb,k));
            x(1:varb,k+1)=c*x(1:varb,k)+d*v(1:varb,k+1);

            x(1:varb,k+1)=feval(@xModifiedPSO,x(1:varb,k+1),lowerbound,upperbound);
            %         fff(k+1)=feval(@func,x(1:varb,k+1))
            fff(k+1)=feval(@TettaHandle,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,x(1:varb,k+1));

            e=e+1;
            p(k)=min(fff);
        end
        xf(:,i)=x(:,end);
        objectf(i)=fff(end);
        min(objectf);
    end
    [ddd,hhh]=min(objectf);
    bestVariables=xf(:,hhh);
    bestFunction=objectf(hhh);
    e;
%     toc;
    set(handles.pw, 'Visible', 'off');
    seconds=toc;

    % Results

    set(handles.not, 'String', bestVariables(1) );
    set(handles.ft, 'String',bestVariables(2)  );
    if bestVariables(3)==1
        set(handles.ct, 'String',  'Total');
    else
        set(handles.ct, 'String',  'Partial');
    end
    set(handles.tc, 'String',bestFunction);
    set(handles.fe, 'String',e );
    set(handles.t, 'String', seconds  );

end



function run_uipanel7_vertical_Callback(hObject, eventdata, handles)


formul=get(handles.formul,'String');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ra=get(handles.ra,'Value');
ga=get(handles.ga,'Value');
pso=get(handles.pso,'Value');


if ra==1

    % parameters
    nov_uipanel9=str2double(get(handles.nov_uipanel9,'String'));
    pp_uipanel9=str2double(get(handles.pp_uipanel9,'String'));
    nop_uipanel9=str2double(get(handles.nop_uipanel9,'String'));
    noi_uipanel9=str2double(get(handles.noi_uipanel9,'String'));
    norp_uipanel9=str2double(get(handles.norp_uipanel9,'String'));
    nord_uipanel9=str2double(get(handles.nord_uipanel9,'String'));
    mint_uipanel9=str2double(get(handles.mint_uipanel9,'String'));
    maxt_uipanel9=str2double(get(handles.maxt_uipanel9,'String'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.pw_vertical, 'Visible', 'on');

    pause(.001)
    %%%%%%%%%%%
    format short
    tic
    %-------------------------------------------------------------------------
    condition(1:nov_uipanel9)={'real'}; % type of variable(s)
    %-------------------------------------------------------------------------
    varb=nov_uipanel9;%   3     ;              % number of variable(s)
    nc=     1     ;              % number of combination between   1  to (varb-2)
    %-------------------------------------------------------------------------
    spp=pp_uipanel9;%    50   ;               % size of primary population (( More than 10 ))
    %-------------------------------------------------------------------------
    np(1:varb)=nop_uipanel9;%   11   ;         % number of parts for real numbers
    npfirst=np;
    secondp =noi_uipanel9;%     1   ;          % number of intervals in " axnew "
    %-------------------------------------------------------------------------
    randomno=  1    ;            % number of steps in Randomizers
    %-------------------------------------------------------------------------
    NOP=norp_uipanel9;%    np(1) ;              % number of parts in Rain Algorithm(Rain Parts)
    N=nord_uipanel9;%      9     ;              % number Of Raindrop
    rp=     0     ;              % number of randoms
    %-------------------------------------------------------------------------
    lowerbound(1:varb)=    mint_uipanel9    ;    % lower bound of variable(s)
    upperbound(1:varb)=    maxt_uipanel9  ;    % upper bound of variable(s)
    %-------------------------------------------------------------------------
    e=1 ;                        % function evaluation
    %-------------------------------------------------------------------------
    for i=1:varb
        realcond=strcmpi(condition(i),'real');
        if realcond == 1
            answer=linspace(lowerbound(i),upperbound(i),np(i));
            ax(i,1:numel(answer))=answer;
        end
        integercond=strcmpi(condition(i),'integer');
        if integercond == 1
            answer=ceil(lowerbound(i)):1:floor(upperbound(i));
            ax(i,1:numel(answer))=answer;
            np(i)=numel(answer);
        end
    end
    ax;
    axfirst=ax;


    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %                  MatrixAlgorithmVersion3 Func : MAVersion3
    [bestfxnew,bestmxnew,bestxnew,ax,np,e]=feval(@MAVersion3FuncOpt2nd,formul,ax,spp,varb,...
        secondp,np,e,nc,randomno,NOP,N,rp,lowerbound,upperbound);
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    bestfxnew';
    bestmxnew;
    bestxnew(:,end)';

    bestfxfirst=bestfxnew(end);
    bestmxfirst=bestmxnew(:,end);
    bestxfirst=bestxnew(:,end)';
    e;


    seconds=toc;

    set(handles.pw_vertical, 'Visible', 'off');

    % Results
    set(handles.x, 'Max',numel(bestxfirst));
    set(handles.x, 'String', bestxfirst);
    set(handles.mv_uipanel13, 'String',  bestfxfirst);
    set(handles.fe_uipanel13, 'String', e  );
    set(handles.t_uipanel13, 'String', seconds  );



elseif ga==1

    % parameters
    nov_uipanel10=str2double(get(handles.nov_uipanel10,'String'));
    nop_uipanel10=str2double(get(handles.nop_uipanel10,'String'));
    nog_uipanel10=str2double(get(handles.nog_uipanel10,'String'));
    nomcg_uipanel10=str2double(get(handles.nomcg_uipanel10,'String'));
    nomcr_uipanel10=str2double(get(handles.nomcr_uipanel10,'String'));
    noec_uipanel10=str2double(get(handles.noec_uipanel10,'String'));
    mint_uipanel10=str2double(get(handles.mint_uipanel10,'String'));
    maxt_uipanel10=str2double(get(handles.maxt_uipanel10,'String'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.pw_vertical, 'Visible', 'on');
    pause(.001)

    tic
    % clc
    % figure(1)
    % clf
    %     clear all
    format short

    %------------------------        parameters        ------------------------
    % befor using this function you must specified your function in fun00.m
    % file in current directory and then set the parameters
    var=nov_uipanel10;%3;            % Number of variables (this item must be equal to the
    %   number of variables that is used in the function in
    %   fun00.m file)
    n=nop_uipanel10;%100;            % Number of population

    m0=nog_uipanel10;%20;            % Number of generations that max value remains constant
    %   (use for termination criteria)
    nmutationG=nomcg_uipanel10;%20;                  %number of mutation children(Gaussian)
    nmutationR=nomcr_uipanel10;%20;                  %number of mutation children(random)
    nelit=noec_uipanel10;%2;                        %number of elitism children
    valuemin=ones(1,var)*   mint_uipanel10 ;      % min possible value of variables
    valuemax=ones(1,var)*  maxt_uipanel10 ;     % max possible value of variables

    %-------------------------------------------------------------------------
    nmutation=nmutationG+nmutationR;
    sigma=(valuemax-valuemin)/10;    %Parameter that related to Gaussian
    %   function and used in mutation step
    max1=zeros(nelit,var);
    parent=zeros(n,var);
    cu=[valuemin(1) valuemax(1) valuemin(2) valuemax(2)];
    for l=1:var
        p(:,l)=valuemin(l)+rand(n,1).*(valuemax(l)-valuemin(l));
    end
    initial=p;
    m=m0;
    maxvalue=ones(m,1)*-1e10;
    maxvalue(m)=-1e5;
    g=0;
    meanvalue(m)=0;
    %-------------   ****    termination criteria   ****-------------
    while abs(maxvalue(m)-maxvalue(m-(m0-1)))>0.001*maxvalue(m) &...
            (abs(maxvalue(m))>1e-10 & abs(maxvalue(m-(m0-1)))>1e-10)...
            & m<10000 & abs(maxvalue(m)-meanvalue(m))>1e-5 | m<20
        sigma=sigma./(1.05);% reducing the sigma value
        %  ------  **** % reducing the number of mutation()random   **** ----
        g=g+1;
        if g>10 & nmutationR>0
            g=0;
            nmutationR=nmutationR-1;
            nmutation=nmutationG+nmutationR;
        end

        %-------------   ****    function evaluation   ****-------------

        for i=1:n
            y(i)=-funcformul(formul,p(i,:));


        end
        s=sort(y);
        maxvalue1(1:nelit)=s(n:-1:n-nelit+1);
        if nelit==0
            maxvalue1(1)=s(n);
            for i=1:n
                if y(i)==maxvalue1(1)
                    max1(1,:)=p(i,:);
                end
            end
        end
        for k=1:nelit
            for i=1:n
                if y(i)==maxvalue1(k)
                    max1(k,:)=p(i,:);
                end
            end
        end


        y=y-min(y)*1.02;
        sumd=y./sum(y);
        meanvalue=y./(sum(y)/n);


        %-------------   ****   Selection: Rolette wheel   ****-------------
        for l=1:n
            sel=rand;
            sumds=0;
            j=1;
            while sumds<sel
                sumds=sumds+sumd(j);
                j=j+1;
            end
            parent(l,:)=p(j-1,:);
        end
        p=zeros(n,var);

        %-------------   ****    regeneration   ****-------------
        for l=1:var


            %-------------   ****    cross-over   ****-------------
            for j=1:ceil((n-nmutation-nelit)/2)
                t=rand*1.5-0.25;
                p(2*j-1,l)=t*parent(2*j-1,l)+(1-t)*parent(2*j,l);
                p(2*j,l)=t*parent(2*j,l)+(1-t)*parent(2*j-1,l);
            end


            %-------------   ****    elitism   ****-------------
            for k=1:nelit
                p((n-nmutation-k+1),l)=max1(k,l);
            end


            %-------------   ****    mutation   ****-------------
            for i=n-nmutation+1:n-nmutationR
                phi=1-2*rand;
                z=erfinv(phi)*(2^0.5);
                p(i,l)=z*sigma(l)+parent(i,l);

            end
            for i=n-nmutationR+1:n
                p(i,1:var)=valuemin(1:var)+rand(1,var).*(valuemax(1:var)-...
                    valuemin(1:var));
            end
            for i=1:n
                for l=1:var
                    if p(i,l)<valuemin(l)
                        p(i,l)=valuemin(l);
                    elseif p(i,l)>valuemax(l)
                        p(i,l)=valuemax(l);
                    end
                end
            end
        end
        p;
        m=m+1;
        max1;
        maxvalue(m)=maxvalue1(1);
        maxvalue00(m-m0)=-maxvalue1(1);
        mean00(m-m0)=sum(s)/n;
        meanvalue(m)=mean00(m-m0);

    end



%     disp('**************************************')
    num_of_fun_evaluation=n*m;
    max_point_GA=max1(1,:);
    maxvalue_GA=maxvalue00(m-m0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    seconds=toc;
    set(handles.pw_vertical, 'Visible', 'off');

    % Results
    set(handles.x, 'Max',numel(max_point_GA));
    set(handles.x, 'String', max_point_GA);
    set(handles.mv_uipanel13, 'String', maxvalue_GA);
    set(handles.fe_uipanel13, 'String', num_of_fun_evaluation  );
    set(handles.t_uipanel13, 'String', seconds  );


elseif pso


    set(handles.pw_vertical, 'Visible', 'on');
    pause(.001)

    % parameters
    nov_uipanel11=str2double(get(handles.nov_uipanel11,'String'));
    pa_uipanel11=(get(handles.pa_uipanel11,'String'));
    pb_uipanel11=(get(handles.pb_uipanel11,'String'));
    pc_uipanel11=str2double(get(handles.pc_uipanel11,'String'));
    pd_uipanel11=str2double(get(handles.pd_uipanel11,'String'));
    ivov_uipanel11=str2double(get(handles.ivov_uipanel11,'String'));
    ivox_uipanel11=str2double(get(handles.ivox_uipanel11,'String'));
    noi_uipanel11=str2double(get(handles.noi_uipanel11,'String'));
    mint_uipanel11=str2double(get(handles.mint_uipanel11,'String'));
    maxt_uipanel11=str2double(get(handles.maxt_uipanel11,'String'));


    parametra=strcmpi(pa_uipanel11,'set');
    parametrb=strcmpi(pb_uipanel11,'set');
    if parametra==1 && parametrb==1
        a=[.9 .7 .9 .1 .1 -.7 ];
        b=[.1 .3  3 .1 2.1 .5 ];
    else
        pa_uipanel11=str2double(get(handles.pa_uipanel11,'String'));
        pb_uipanel11=str2double(get(handles.pb_uipanel11,'String'));
        a=pa_uipanel11;%[.9 .7 .9 .1 .1 -.7 ];
        b=pb_uipanel11;%[.1 .3  3 .1 2.1 .5 ];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    format short
    tic
    e=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varb=nov_uipanel11;%3;

    c=pc_uipanel11;%;
    d=pd_uipanel11;%1;

    lowerbound=mint_uipanel11;%   4   ;    % lower bound of variable(s)
    upperbound= maxt_uipanel11;%  100  ;    % upper bound of variable(s)

    for i=1:numel(a)

        v(1:varb,1)=ivov_uipanel11;%-0.1;
        x(1:varb,1)=ivox_uipanel11;%10;
        p(1)=feval(@funcformul,formul,x(1:varb,1));

        v(1:varb,2)=a(i)*v(1:varb,1)+b(i)*(p(1)-x(1:varb,1));
        x(1:varb,2)=c*x(1:varb,1)+d*v(1:varb,1);

        k=1;
        for j=1:noi_uipanel11;%500
            k=k+1;
            v(1:varb,k+1)=a(i)*v(1:varb,k)+b(i)*rand(varb,1).*(p(k-1)-x(1:varb,k));
            x(1:varb,k+1)=c*x(1:varb,k)+d*v(1:varb,k+1);

            for kk=1:numel(x)
                if x(kk)<lowerbound || x(kk)>upperbound
                    x(kk)=round(lowerbound+rand*(upperbound-lowerbound));
                end
            end
            %                         x(1:varb,k+1)=feval(@xModifiedPSO,x(1:varb,k+1),lowerbound,upperbound);
            fff(k+1)=feval(@funcformul,formul,x(1:varb,k+1));
            %             fff(k+1)=feval(@TettaHandle,x(1:varb,k+1));

            e=e+1;
            p(k)=min(fff);
        end
        xf(:,i)=x(:,end);
        objectf(i)=fff(end);
        min(objectf);
    end
    [ddd,hhh]=min(objectf);
    bestVariables=xf(:,hhh);
    bestFunction=objectf(hhh);
    e;
    %toc;
    set(handles.pw_vertical, 'Visible', 'off');
    seconds=toc;

    % Results
    set(handles.x, 'Max',numel(bestVariables));
    set(handles.x, 'String', bestVariables);
    set(handles.mv_uipanel13, 'String',bestFunction);
    set(handles.fe_uipanel13, 'String',e  );
    set(handles.t_uipanel13, 'String', seconds  );

end


function reset_uipanel7_Callback(hObject, eventdata, handles)

set(handles.tc, 'String','');
set(handles.not, 'String','');
set(handles.ft, 'String','');
set(handles.fe, 'String','');
set(handles.t, 'String','');
set(handles.ct, 'String','');
set(handles.pw, 'Visible', 'off');

set(handles.fr_uipanel23, 'String','');
set(handles.d_uipanel23, 'String','');
set(handles.r_uipanel23, 'String','');
set(handles.ftemp_uipanel23, 'String','');
set(handles.bldlupper, 'String','');
set(handles.bhdhlower, 'String','');

set(handles.fs_uipanel23, 'Value', 1 );
set(handles.noc_uipanel23, 'Value', 1 );
set(handles.uipanel25, 'Visible', 'off');

set(handles.mfoc11 , 'String', '' );
set(handles.mfoc22 , 'String',  '' );
set(handles.mfoc33 , 'String',  '' );
set(handles.mfoc44, 'String', '');
set(handles.mfoc55, 'String', '');
set(handles.mfoc66, 'String', '');
set(handles.mfoc77, 'String', '');
set(handles.mfoc88, 'String', '');
set(handles.mfoc99, 'String', '');
set(handles.mfoc1010, 'String', '');

function reset_uipanel7_vertical_Callback(hObject, eventdata, handles)

set(handles.x, 'String','');
set(handles.mv_uipanel13, 'String','');
set(handles.fe_uipanel13, 'String','');
set(handles.t_uipanel13, 'String','');
set(handles.pw_vertical, 'Visible', 'off');



%%%%%%%%%%%%%%%    dp_uipanel7 Parameters Button_uipanel7

function dp_uipanel7_Callback(hObject, eventdata, handles)

ra=get(handles.ra,'Value');
ga=get(handles.ga,'Value');
pso=get(handles.pso,'Value');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ra==1
    % Rain (defalt parameters)
    set(handles.nov_uipanel9,'String','3');
    set(handles.pp_uipanel9,'String','50');
    set(handles.nop_uipanel9,'String','11');
    set(handles.noi_uipanel9,'String','1');
    set(handles.norp_uipanel9,'String','11');
    set(handles.nord_uipanel9,'String','9');
    set(handles.mint_uipanel9,'String','4');
    set(handles.maxt_uipanel9,'String','100');

    set(handles.uipanel10, 'Visible', 'off');
    set(handles.uipanel11, 'Visible', 'off');

    set(handles.uipanel9, 'Visible', 'on');
    set(handles.uipanel9, 'Position', [31.2,29,39.6,25.846]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ga==1
    % GA (defalt parameters)
    set(handles.nov_uipanel10,'String','3');
    set(handles.nop_uipanel10,'String','100');
    set(handles.nog_uipanel10,'String','20');
    set(handles.nomcg_uipanel10,'String','20');
    set(handles.nomcr_uipanel10,'String','20');
    set(handles.noec_uipanel10,'String','2');
    set(handles.mint_uipanel10,'String','4');
    set(handles.maxt_uipanel10,'String','100');

    set(handles.uipanel9, 'Visible', 'off');
    set(handles.uipanel11, 'Visible', 'off');

    set(handles.uipanel10, 'Visible', 'on');
    set(handles.uipanel10, 'Position', [31.2,29,39.6,25.846]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pso==1
    % PSO (defalt parameters)
    set(handles.nov_uipanel11,'String','3');
    set(handles.pa_uipanel11,'String','set');
    set(handles.pb_uipanel11,'String','set');
    set(handles.pc_uipanel11,'String','1');
    set(handles.pd_uipanel11,'String','1');
    set(handles.ivov_uipanel11,'String','-0.1');
    set(handles.ivox_uipanel11,'String','1');
    set(handles.noi_uipanel11,'String','500');
    set(handles.mint_uipanel11,'String','4');
    set(handles.maxt_uipanel11,'String','100');


    set(handles.uipanel9, 'Visible', 'off');
    set(handles.uipanel10, 'Visible', 'off');
    %
    set(handles.uipanel11, 'Visible', 'on');
    set(handles.uipanel11, 'Position', [31.2,29,39.6,25.846]);
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------































function tc_Callback(hObject, eventdata, handles)

function tc_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function not_Callback(hObject, eventdata, handles)

function not_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ft_Callback(hObject, eventdata, handles)

function ft_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fe_Callback(hObject, eventdata, handles)

function fe_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function listbox1_Callback(hObject, eventdata, handles)

function listbox1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider1_Callback(hObject, eventdata, handles)

function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popupmenu1_Callback(hObject, eventdata, handles)

function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function togglebutton2_Callback(hObject, eventdata, handles)

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_Callback(hObject, eventdata, handles)

function t_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function password_Callback(hObject, eventdata, handles)

function password_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in popupmenu1.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in back_uipanel23.





function formul_Callback(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of formul as text
%        str2double(get(hObject,'String')) returns contents of formul as a double


% --- Executes during object creation, after setting all properties.
function formul_CreateFcn(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of formul as text
%        str2double(get(hObject,'String')) returns contents of formul as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok.




% --- Executes on button press in stp_uipanel7.
function stp_uipanel7_Callback(hObject, eventdata, handles)
% hObject    handle to stp_uipanel7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)








function nov_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to nov_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nov_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nov_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function nov_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nov_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pp_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to pp_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pp_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of pp_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function pp_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pp_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nop_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to nop_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nop_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of nop_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function nop_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nop_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noi_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to noi_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noi_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of noi_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function noi_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noi_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function norp_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to norp_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of norp_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of norp_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function norp_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to norp_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nord_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to nord_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nord_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of nord_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function nord_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nord_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mint_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to mint_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mint_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of mint_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function mint_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mint_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxt_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxt_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of maxt_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function maxt_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to nov_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nov_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nov_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nov_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nop_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to nop_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nop_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nop_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function nop_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nop_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nog_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to nog_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nog_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nog_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function nog_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nog_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nomcg_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to nomcg_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nomcg_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nomcg_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function nomcg_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nomcg_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nomcr_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to nomcr_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nomcr_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of nomcr_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function nomcr_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nomcr_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noec_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to noec_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noec_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of noec_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function noec_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noec_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mint_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to mint_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mint_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of mint_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function mint_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mint_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxt_uipanel10_Callback(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxt_uipanel10 as text
%        str2double(get(hObject,'String')) returns contents of maxt_uipanel10 as a double


% --- Executes during object creation, after setting all properties.
function maxt_uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nov_uipanel9_Callback(hObject, eventdata, handles)
% hObject    handle to nov_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nov_uipanel9 as text
%        str2double(get(hObject,'String')) returns contents of nov_uipanel9 as a double


% --- Executes during object creation, after setting all properties.
function nov_uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nov_uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function nov_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to nov_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nov_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of nov_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function nov_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nov_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pa_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to pa_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pa_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of pa_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function pa_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pa_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pb_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to pb_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pb_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of pb_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function pb_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pb_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pc_uipanel12_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to pc_uipanel12_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc_uipanel12_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of pc_uipanel12_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function pc_uipanel12_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc_uipanel12_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to pd_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of pd_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function pd_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ivov_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to ivov_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ivov_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of ivov_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function ivov_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ivov_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ivox_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to ivox_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ivox_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of ivox_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function ivox_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ivox_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noi_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to noi_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noi_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of noi_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function noi_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noi_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function mint_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to mint_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mint_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of mint_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function mint_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mint_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxt_uipanel11_Callback(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxt_uipanel11 as text
%        str2double(get(hObject,'String')) returns contents of maxt_uipanel11 as a double


% --- Executes during object creation, after setting all properties.
function maxt_uipanel11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxt_uipanel11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ct_Callback(hObject, eventdata, handles)
% hObject    handle to ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ct as text
%        str2double(get(hObject,'String')) returns contents of ct as a double


% --- Executes during object creation, after setting all properties.
function ct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ft as text
%        str2double(get(hObject,'String')) returns contents of ft as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ct as text
%        str2double(get(hObject,'String')) returns contents of ct as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tc as text
%        str2double(get(hObject,'String')) returns contents of tc as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fe as text
%        str2double(get(hObject,'String')) returns contents of fe as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t as text
%        str2double(get(hObject,'String')) returns contents of t as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function list_Callback(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of list as text
%        str2double(get(hObject,'String')) returns contents of list as a double


% --- Executes during object creation, after setting all properties.
function list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit55_Callback(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit55 as text
%        str2double(get(hObject,'String')) returns contents of edit55 as a double


% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit59_Callback(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59 as text
%        str2double(get(hObject,'String')) returns contents of edit59 as a double


% --- Executes during object creation, after setting all properties.
function edit59_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit60_Callback(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60 as text
%        str2double(get(hObject,'String')) returns contents of edit60 as a double


% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit61_Callback(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit61 as text
%        str2double(get(hObject,'String')) returns contents of edit61 as a double


% --- Executes during object creation, after setting all properties.
function edit61_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit62_Callback(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit62 as text
%        str2double(get(hObject,'String')) returns contents of edit62 as a double


% --- Executes during object creation, after setting all properties.
function edit62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit63_Callback(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit63 as text
%        str2double(get(hObject,'String')) returns contents of edit63 as a double


% --- Executes during object creation, after setting all properties.
function edit63_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit64_Callback(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64 as text
%        str2double(get(hObject,'String')) returns contents of edit64 as a double


% --- Executes during object creation, after setting all properties.
function edit64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit65 as text
%        str2double(get(hObject,'String')) returns contents of edit65 as a double


% --- Executes during object creation, after setting all properties.
function edit65_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in reset_uipanel7_vertical.


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit72_Callback(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of formul as text
%        str2double(get(hObject,'String')) returns contents of formul as a double


% --- Executes during object creation, after setting all properties.
function edit72_CreateFcn(hObject, eventdata, handles)
% hObject    handle to formul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


















function x_Callback(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double


% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function ok_Callback(hObject, eventdata, handles)

% formul=get(handles.formul,'String');
% formul;
% var=[8 7];
% answer=feval(@testformul,var,formul);

% --- Executes on button press in hmd.
% function hmd_Callback(hObject, eventdata, handles)
% x=[1  2  3  4 6 6 6 6 6 7 7 7 7 8 8 8 8 8 ];
% set(handles.list,'String',x);
%












function mv_uipanel13_Callback(hObject, eventdata, handles)
% hObject    handle to mv_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mv_uipanel13 as text
%        str2double(get(hObject,'String')) returns contents of mv_uipanel13 as a double


% --- Executes during object creation, after setting all properties.
function mv_uipanel13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mv_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fe_uipanel13_Callback(hObject, eventdata, handles)
% hObject    handle to fe_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fe_uipanel13 as text
%        str2double(get(hObject,'String')) returns contents of fe_uipanel13 as a double


% --- Executes during object creation, after setting all properties.
function fe_uipanel13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fe_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_uipanel13_Callback(hObject, eventdata, handles)
% hObject    handle to t_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_uipanel13 as text
%        str2double(get(hObject,'String')) returns contents of t_uipanel13 as a double


% --- Executes during object creation, after setting all properties.
function t_uipanel13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next_uipanel6.



% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit83_Callback(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit83 as text
%        str2double(get(hObject,'String')) returns contents of edit83 as a double


% --- Executes during object creation, after setting all properties.
function edit83_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit84_Callback(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit84 as text
%        str2double(get(hObject,'String')) returns contents of edit84 as a double


% --- Executes during object creation, after setting all properties.
function edit84_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit85_Callback(hObject, eventdata, handles)
% hObject    handle to edit85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit85 as text
%        str2double(get(hObject,'String')) returns contents of edit85 as a double


% --- Executes during object creation, after setting all properties.
function edit85_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit86_Callback(hObject, eventdata, handles)
% hObject    handle to edit86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit86 as text
%        str2double(get(hObject,'String')) returns contents of edit86 as a double


% --- Executes during object creation, after setting all properties.
function edit86_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit87_Callback(hObject, eventdata, handles)

function edit87_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit88_Callback(hObject, eventdata, handles)

function edit88_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pushbutton39_Callback(hObject, eventdata, handles)

function pushbutton42_Callback(hObject, eventdata, handles)

function run_uipanel12_Callback(hObject, eventdata, handles)

set(handles.pw_uipanel12,'Visible','on');
pause(0.001)
popup_sel_index = get(handles.fs, 'Value');
switch popup_sel_index
    case 2
        cp=str2double(get(handles.cp,'String'));
        bp=str2double(get(handles.bp,'String'));
        tf=str2double(get(handles.tf,'String'));
        lh=str2double(get(handles.lh,'String'));
        q=1+cp*(bp-tf)/lh;
    case 3
        q=1;
    case 4
        q=str2double(get(handles.mfl,'String'));

    case 5
        q=0;
    case 6
        cp=str2double(get(handles.cp,'String'));
        bp=str2double(get(handles.bp,'String'));
        tf=str2double(get(handles.tf,'String'));
        lh=str2double(get(handles.lh,'String'));
        q=1+cp*(bp-tf)/lh;
end

%read from input
f=str2double(get(handles.f,'String'));
zf=str2double(get(handles.zf,'String'));
xd=str2double(get(handles.xd,'String'));
xw=str2double(get(handles.xw,'String'));
alpha=str2double(get(handles.alpha,'String'));
r=str2double(get(handles.r,'String'));

%calculation of w , d , rm, Nm,
w=f*(zf-xd)/(xw-xd);
d=f-w;

Nm=log10(xd*(1-xw)/(xw*(1-xd)))/log10(alpha);

%calculaton of N & x , y in each tray
l=r*d;
v=(r+1)*d;
lb=l+f*q;
vb=v+f*(q-1);
xe=inline('y/(alpha-y*(alpha-1))','alpha','y');
yu=inline('r/(r+1)*x+xd/(r+1)','r','xd','x');
yb=inline('lb/(lb-w)*x-w*xw/(lb-w)','lb','w','xw','x');
yf=inline('q/(q-1)*x-zf/(q-1)','q','zf');
% equlibrium condensor

xcross=inline('(xd*q-xd+zf*r+zf)/(r+q)','q','r','zf','xd');
xcross=xcross(q,r,zf,xd);
i=1;
tc=get(handles.tc_uipanel12,'Value');
if tc==1
    xx=xd;
else
    xx=xe(alpha,xd);
end
y(i)=yu(r,xd,xx);
x(i)=xe(alpha,y(i));
i=i+1;


while xx>xcross
    y(i)=yu(r,xd,x(i-1));
    x(i)=xe(alpha,y(i));
    xx=x(i);
    i=i+1;


end

while xx>xw
    y(i)=yb(lb,w,xw,x(i-1));
    x(i)=xe(alpha,y(i));
    xx=x(i);
    i=i+1;
end
t=(x(i-2)-xw)/(x(i-2)-x(i-1));
%nu of tray
for j=1:(length(x)-1)
    mmm(j)=x(j);
    nnn(j)=y(j);
end
if tc==1
    N=length(mmm)-1+t;
else
    N=length(nnn)-1+t;
end
set(handles.nort, 'String', round(N));
W=f/(1+((xw-zf)/(zf-xd)));
D=f-W;
set(handles.D, 'String', D);
set(handles.W, 'String', W);

set(handles.pw_uipanel12,'Visible','off');



function f_Callback(hObject, eventdata, handles)

function f_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xd_Callback(hObject, eventdata, handles)

function xd_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xw_Callback(hObject, eventdata, handles)

function xw_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zf_Callback(hObject, eventdata, handles)

function zf_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)

function alpha_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_Callback(hObject, eventdata, handles)

function r_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refluxlimit.
function refluxlimit_Callback(hObject, eventdata, handles)

set(handles.pw_uipanel12,'Visible','on');
pause(0.001)
popup_sel_index = get(handles.fs, 'Value');
switch popup_sel_index
    case 2
        cp=str2double(get(handles.cp,'String'));
        bp=str2double(get(handles.bp,'String'));
        tf=str2double(get(handles.tf,'String'));
        lh=str2double(get(handles.lh,'String'));
        q=1+cp*(bp-tf)/lh;
    case 3
        q=1;
    case 4
        q=str2double(get(handles.mfl,'String'));

    case 5
        q=0;
    case 6
        cp=str2double(get(handles.cp,'String'));
        bp=str2double(get(handles.bp,'String'));
        tf=str2double(get(handles.tf,'String'));
        lh=str2double(get(handles.lh,'String'));
        q=1+cp*(bp-tf)/lh;
end
zf=str2double(get(handles.zf,'String'));
xd=str2double(get(handles.xd,'String'));
xw=str2double(get(handles.xw,'String'));
alpha=str2double(get(handles.alpha,'String'));
r=str2double(get(handles.r,'String'));

rm=0;
f1=1;f2=0;
while abs(f1-f2)>0.05
    rm=rm+.05;
    f1=(rm*zf+q*xd)/(rm*(1-zf)+q*(1-xd));
    f2=(alpha*(xd*(q-1)+zf*(rm+1)))/((rm+1)*(1-zf)+(q-1)*(1-xd));
end
set(handles.minreflux, 'String', 1.2*rm);
set(handles.maxreflux, 'String', 1.5*rm);

set(handles.pw_uipanel12,'Visible','off');




function dp_uipanel12_Callback(hObject, eventdata, handles)

set(handles.f , 'String',  216.8  );
set(handles.xd , 'String',  0.915  );
set(handles.xw , 'String',  0.00565  );
set(handles.zf , 'String',  0.36  );
set(handles.alpha , 'String',  4.16  );
set(handles.r , 'String', 2.25   );
set(handles.mfl , 'String', 0.36  );
set(handles.tc_uipanel12 , 'Value', 1  );
set(handles.fs , 'Value', 4  );
set(handles.pw_uipanel12, 'Visible', 'off');

set(handles.feedpanel, 'Visible', 'off');
set(handles.molefractionpanel, 'Visible', 'on');



function reset_uipanel12_Callback(hObject, eventdata, handles)

set(handles.f, 'String','');
set(handles.xd, 'String','');
set(handles.xw, 'String','');
set(handles.zf, 'String','');
set(handles.alpha, 'String','');
set(handles.r, 'String','');

set(handles.cp, 'String','');
set(handles.tf, 'String','');
set(handles.bp, 'String','');
set(handles.lh, 'String','');
set(handles.mfl, 'String','');
set(handles.nort, 'String','');
set(handles.minreflux, 'String','');
set(handles.maxreflux, 'String','');
set(handles.D, 'String','');
set(handles.W, 'String','');

set(handles.pw_uipanel12, 'Visible', 'off');


function pushbutton47_Callback(hObject, eventdata, handles)




function cp_Callback(hObject, eventdata, handles)

function cp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bp_Callback(hObject, eventdata, handles)

function bp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tf_Callback(hObject, eventdata, handles)

function tf_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lh_Callback(hObject, eventdata, handles)

function lh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_Callback(hObject, eventdata, handles)

set(handles.feedpanel, 'Visible', 'off');
set(handles.molefractionpanel, 'Visible', 'off');



popup_sel_index = get(handles.fs, 'Value');
switch popup_sel_index
    case 1
        
        set(handles.feedpanel, 'Visible', 'off');
        set(handles.molefractionpanel, 'Visible', 'off');
    
    case 2
        
        set(handles.feedpanel, 'Visible', 'on');
        set(handles.molefractionpanel, 'Visible', 'off');
    
    case 3
        
        set(handles.feedpanel, 'Visible', 'off');
        set(handles.molefractionpanel, 'Visible', 'off');
    
    case 4
        
        set(handles.molefractionpanel, 'Visible', 'on');
        set(handles.feedpanel, 'Visible', 'off');
    
    case 5
        
        set(handles.feedpanel, 'Visible', 'off');
        set(handles.molefractionpanel, 'Visible', 'off');
   
    case 6
        
        set(handles.feedpanel, 'Visible', 'on');
        set(handles.molefractionpanel, 'Visible', 'off');
        
end


function fs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfl_Callback(hObject, eventdata, handles)

function mfl_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nort_Callback(hObject, eventdata, handles)

function nort_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pc_uipanel11_Callback(hObject, eventdata, handles)


function pc_uipanel11_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit110_Callback(hObject, eventdata, handles)
% hObject    handle to edit110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit110 as text
%        str2double(get(hObject,'String')) returns contents of edit110 as a double


% --- Executes during object creation, after setting all properties.
function edit110_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit110 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit111_Callback(hObject, eventdata, handles)
% hObject    handle to edit111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit111 as text
%        str2double(get(hObject,'String')) returns contents of edit111 as a double


% --- Executes during object creation, after setting all properties.
function edit111_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit111 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit112_Callback(hObject, eventdata, handles)
% hObject    handle to edit112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit112 as text
%        str2double(get(hObject,'String')) returns contents of edit112 as a double


% --- Executes during object creation, after setting all properties.
function edit112_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit112 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit113_Callback(hObject, eventdata, handles)
% hObject    handle to edit113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit113 as text
%        str2double(get(hObject,'String')) returns contents of edit113 as a double


% --- Executes during object creation, after setting all properties.
function edit113_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit113 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fr_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to fr_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of fr_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function fr_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to d_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of d_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function d_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to r_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of r_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function r_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit107_Callback(hObject, eventdata, handles)
% hObject    handle to edit107 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit107 as text
%        str2double(get(hObject,'String')) returns contents of edit107 as a double


% --- Executes during object creation, after setting all properties.
function edit107_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit107 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ft_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to ft_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ft_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of ft_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function ft_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ft_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit109_Callback(hObject, eventdata, handles)
% hObject    handle to edit109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit109 as text
%        str2double(get(hObject,'String')) returns contents of edit109 as a double


% --- Executes during object creation, after setting all properties.
function edit109_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit109 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton54.
function pushbutton54_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset_uipanel19.
function reset_uipanel19_Callback(hObject, eventdata, handles)

set(handles.fr_uipanel19, 'String','');
set(handles.xd_uipanel19, 'String','');
set(handles.xw_uipanel19, 'String','');
set(handles.d_uipanel19, 'String','');
set(handles.r_uipanel19, 'String','');
set(handles.ft_uipanel19, 'String','');
set(handles.ftemp_uipanel19, 'String','');
set(handles.not_uipanel19, 'String','');
set(handles.noc, 'Value', 1 );
set(handles.fs_uipanel19, 'Value',1);
set(handles.qc, 'String','');
set(handles.qr, 'String','');

set(handles.uipanel21, 'Visible', 'off');

set(handles.mfoc1, 'String','')
set(handles.mfoc2, 'String','')
set(handles.mfoc3, 'String','')
set(handles.mfoc4, 'String','')
set(handles.mfoc5, 'String','')
set(handles.mfoc6, 'String','')
set(handles.mfoc7, 'String','')
set(handles.mfoc8, 'String','')
set(handles.mfoc9, 'String','')
set(handles.mfoc10, 'String','')


set(handles.pw_uipanel19, 'Visible', 'off');

% --- Executes on button press in run_uipanel19.
function run_uipanel19_Callback(hObject, eventdata, handles)

set(handles.pw_uipanel19,'Visible','on');
pause(.001)

fr_uipanel19=str2double(get(handles.fr_uipanel19,'String'));
d_uipanel19=str2double(get(handles.d_uipanel19,'String'));
r_uipanel19=str2double(get(handles.r_uipanel19,'String'));
ftemp_uipanel19=str2double(get(handles.ftemp_uipanel19,'String'));
not_uipanel19=str2double(get(handles.not_uipanel19,'String'));
ft_uipanel19=str2double(get(handles.ft_uipanel19,'String'));

mfoc(1)=str2double(get(handles.mfoc1, 'String'));
mfoc(2)=str2double(get(handles.mfoc2, 'String'));
mfoc(3)=str2double(get(handles.mfoc3, 'String'));
mfoc(4)=str2double(get(handles.mfoc4, 'String'));
mfoc(5)=str2double(get(handles.mfoc5, 'String'));
mfoc(6)=str2double(get(handles.mfoc6, 'String'));
mfoc(7)=str2double(get(handles.mfoc7, 'String'));
mfoc(8)=str2double(get(handles.mfoc8, 'String'));
mfoc(9)=str2double(get(handles.mfoc9, 'String'));
mfoc(10)=str2double(get(handles.mfoc10, 'String'));


noc=get(handles.noc, 'Value');
for i=1:noc
    X(i)=mfoc(i);
end

tc_uipanel19=get(handles.tc_uipanel19,'Value');
pc_uipanel19=get(handles.pc_uipanel19,'Value');
fs_uipanel19=get(handles.fs_uipanel19,'Value');

%-----------------------------------------------------------------
format short
tic
%%% Specifications
% bldlupper=   5e-4   ;
% bhdhlower=   3.1   ;

X;           % mole fraction
base= 1 ;                         % base component for calculation of alfa
% base : between { 1,2, ... , c }
F= fr_uipanel19 ;                          % [  lbmol/h  ]
D= d_uipanel19 ;                           % [  lbmol/h  ]
B= F-D ;                          % [  lbmol/h  ]
N= not_uipanel19 ;                            % number of stages
f= ft_uipanel19  ;                            % feed stage
c=numel(X);                       % number of component

%%%  Comment one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tc_uipanel19==1
    condensor={'Total-Condensor'};
elseif pc_uipanel19==1
    condensor={'Partial-Condensor'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fs_uipanel19==2 || fs_uipanel19==3
    feedstate={'Buuble-Point or Subcooled Liquid'};
elseif fs_uipanel19==4 || fs_uipanel19==5
    feedstate={'Dew-Pointe or Super heated Vapor'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initiation

TF= ftemp_uipanel19 ;                         % feed temperature [  "R   ]
T(1:N)= ftemp_uipanel19 ;      % temperature      [  "R   ]

% Assume Vj  ( j= 1, ... ,N )

Vinitial= D*(r_uipanel19+1);

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

%%% Enthalpy functions   %%%   Change it to function @
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

%%% While Loop
for loop=1: 15
    %%   K  values
    %%%             Ki = Ci exp ( -Ei / T )     ,  [T]: "R

    C=[ 4e3 , 8e3 , 12e3  ];
    E=[ 4.6447e3 , 4.6447e3 , 4.6447e3 ];
    for i=1:c
        for j=1:N
            k(j,i)=C(i)*exp(-E(i)/T(j));
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
    end

    SUMalfaXb=sum(alfa(N,:)*Xb(:));
    K(N,base)=1/SUMalfaXb;
    T(N)=E(base)/(log(C(base)/K(N,base)));

    temperature=T';

    %%  Enthalpies   &   Utilities Duty

    for j=1:N
        %     h(j,1)=h1(T(j));
        %     h(j,2)=h2(T(j));
        %     h(j,3)=h3(T(j));
        [h(j,1)]=feval(@h1,T(j));
        [h(j,2)]=feval(@h2,T(j));
        [h(j,3)]=feval(@h3,T(j));
        %     H(j,1)=H1(T(j));
        %     H(j,2)=H2(T(j));
        %     H(j,3)=H3(T(j));
        [H(j,1)]=feval(@HH1,T(j));
        [H(j,2)]=feval(@HH2,T(j));
        [H(j,3)]=feval(@HH3,T(j));
    end
    % Vapor Feed
    %     HF(1)=H1(TF);
    %     HF(2)=H2(TF);
    %     HF(3)=H3(TF);
    [HF(1)]=feval(@HH1,TF);
    [HF(2)]=feval(@HH2,TF);
    [HF(3)]=feval(@HH3,TF);
    % Liquid Feed
    %     hF(1)=h1(TF);
    %     hF(2)=h2(TF);
    %     hF(3)=h3(TF);
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
    % D=sum(d(:))
    % B=sum(b(:))


    %% End While Loop
end
%%%

teta;
b=b';
d=d';

xd=d/D;
xw=b/B;

set(handles.xd_uipanel19, 'String',  xd);
set(handles.xw_uipanel19, 'String', xw  );
set(handles.qr, 'String', Qr  );
set(handles.qc, 'String', Qc);

set(handles.pw_uipanel19,'Visible','off');

%-----------------------------------------------------------------





function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit103_Callback(hObject, eventdata, handles)
% hObject    handle to edit103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit103 as text
%        str2double(get(hObject,'String')) returns contents of edit103 as a double


% --- Executes during object creation, after setting all properties.
function edit103_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit103 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit102_Callback(hObject, eventdata, handles)
% hObject    handle to edit102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit102 as text
%        str2double(get(hObject,'String')) returns contents of edit102 as a double


% --- Executes during object creation, after setting all properties.
function edit102_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit102 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in back_uipanel19.



% --- Executes on button press in dp_uipanel19.
function dp_uipanel19_Callback(hObject, eventdata, handles)

set(handles.fr_uipanel19 , 'String',  100  );
set(handles.d_uipanel19 , 'String',  50  );
set(handles.r_uipanel19 , 'String',  1  );
set(handles.ftemp_uipanel19 , 'String',  400  );
set(handles.not_uipanel19 , 'String',  26 );
set(handles.ft_uipanel19 , 'String', 12   );
set(handles.noc , 'Value',  3  );

set(handles.uipanel21, 'Visible', 'on');
set(handles.mfoc1 , 'String', 0.33  );
set(handles.mfoc2 , 'String', 0.33  );
set(handles.mfoc3 , 'String', 0.33  );
set(handles.mfoc1, 'Visible', 'on');
set(handles.mfoc2, 'Visible', 'on');
set(handles.mfoc3, 'Visible', 'on');
set(handles.mfoctext1, 'Visible', 'on');
set(handles.mfoctext2, 'Visible', 'on');
set(handles.mfoctext3, 'Visible', 'on');

set(handles.mfoc4, 'Visible', 'off');
set(handles.mfoc5, 'Visible', 'off');
set(handles.mfoc6, 'Visible', 'off');
set(handles.mfoc7, 'Visible', 'off');
set(handles.mfoc8, 'Visible', 'off');
set(handles.mfoc9, 'Visible', 'off');
set(handles.mfoc10, 'Visible', 'off');
set(handles.mfoctext4, 'Visible', 'off');
set(handles.mfoctext5, 'Visible', 'off');
set(handles.mfoctext6, 'Visible', 'off');
set(handles.mfoctext7, 'Visible', 'off');
set(handles.mfoctext8, 'Visible', 'off');
set(handles.mfoctext9, 'Visible', 'off');
set(handles.mfoctext10, 'Visible', 'off');

set(handles.tc_uipanel19 , 'Value', 1  );
set(handles.fs_uipanel19 , 'Value', 4  );

set(handles.pw_uipanel19, 'Visible', 'off');




function D_Callback(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D as text
%        str2double(get(hObject,'String')) returns contents of D as a double


% --- Executes during object creation, after setting all properties.
function D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function W_Callback(hObject, eventdata, handles)
% hObject    handle to W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W as text
%        str2double(get(hObject,'String')) returns contents of W as a double


% --- Executes during object creation, after setting all properties.
function W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in noc.
function noc_Callback(hObject, eventdata, handles)


set(handles.uipanel21, 'Visible', 'on');

popup_sel_index = get(handles.noc, 'Value');
switch popup_sel_index

    case 1
        set(handles.mfoc1, 'Visible', 'off');
        set(handles.mfoc2, 'Visible', 'off');
        set(handles.mfoc3, 'Visible', 'off');
        set(handles.mfoc4, 'Visible', 'off');
        set(handles.mfoc5, 'Visible', 'off');
        set(handles.mfoc6, 'Visible', 'off');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
        set(handles.mfoctext1, 'Visible', 'off');
        set(handles.mfoctext2, 'Visible', 'off');
        set(handles.mfoctext3, 'Visible', 'off');
        set(handles.mfoctext4, 'Visible', 'off');
        set(handles.mfoctext5, 'Visible', 'off');
        set(handles.mfoctext6, 'Visible', 'off');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 2
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'off');
        set(handles.mfoc4, 'Visible', 'off');
        set(handles.mfoc5, 'Visible', 'off');
        set(handles.mfoc6, 'Visible', 'off');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
        set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'off');
        set(handles.mfoctext4, 'Visible', 'off');
        set(handles.mfoctext5, 'Visible', 'off');
        set(handles.mfoctext6, 'Visible', 'off');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 3
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'off');
        set(handles.mfoc5, 'Visible', 'off');
        set(handles.mfoc6, 'Visible', 'off');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
         set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'off');
        set(handles.mfoctext5, 'Visible', 'off');
        set(handles.mfoctext6, 'Visible', 'off');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 4
         set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'off');
        set(handles.mfoc6, 'Visible', 'off');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
        set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'off');
        set(handles.mfoctext6, 'Visible', 'off');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 5
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'off');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
        set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'off');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 6
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'on');
        set(handles.mfoc7, 'Visible', 'off');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
         set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'on');
        set(handles.mfoctext7, 'Visible', 'off');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 7
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'on');
        set(handles.mfoc7, 'Visible', 'on');
        set(handles.mfoc8, 'Visible', 'off');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
        set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'on');
        set(handles.mfoctext7, 'Visible', 'on');
        set(handles.mfoctext8, 'Visible', 'off');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 8
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'on');
        set(handles.mfoc7, 'Visible', 'on');
        set(handles.mfoc8, 'Visible', 'on');
        set(handles.mfoc9, 'Visible', 'off');
        set(handles.mfoc10, 'Visible', 'off');
         set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'on');
        set(handles.mfoctext7, 'Visible', 'on');
        set(handles.mfoctext8, 'Visible', 'on');
        set(handles.mfoctext9, 'Visible', 'off');
        set(handles.mfoctext10, 'Visible', 'off');
    case 9
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'on');
        set(handles.mfoc7, 'Visible', 'on');
        set(handles.mfoc8, 'Visible', 'on');
        set(handles.mfoc9, 'Visible', 'on');
        set(handles.mfoc10, 'Visible', 'off');
         set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'on');
        set(handles.mfoctext7, 'Visible', 'on');
        set(handles.mfoctext8, 'Visible', 'on');
        set(handles.mfoctext9, 'Visible', 'on');
        set(handles.mfoctext10, 'Visible', 'off');
    case 10
        set(handles.mfoc1, 'Visible', 'on');
        set(handles.mfoc2, 'Visible', 'on');
        set(handles.mfoc3, 'Visible', 'on');
        set(handles.mfoc4, 'Visible', 'on');
        set(handles.mfoc5, 'Visible', 'on');
        set(handles.mfoc6, 'Visible', 'on');
        set(handles.mfoc7, 'Visible', 'on');
        set(handles.mfoc8, 'Visible', 'on');
        set(handles.mfoc9, 'Visible', 'on');
        set(handles.mfoc10, 'Visible', 'on');
        set(handles.mfoctext1, 'Visible', 'on');
        set(handles.mfoctext2, 'Visible', 'on');
        set(handles.mfoctext3, 'Visible', 'on');
        set(handles.mfoctext4, 'Visible', 'on');
        set(handles.mfoctext5, 'Visible', 'on');
        set(handles.mfoctext6, 'Visible', 'on');
        set(handles.mfoctext7, 'Visible', 'on');
        set(handles.mfoctext8, 'Visible', 'on');
        set(handles.mfoctext9, 'Visible', 'on');
        set(handles.mfoctext10, 'Visible', 'on');

end


% --- Executes during object creation, after setting all properties.
function noc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextcomponent.
function nextcomponent_Callback(hObject, eventdata, handles)
popup_sel_index = get(handles.noc, 'Value');
switch popup_sel_index

    case 2
        c=2;
    case 3
        c=3;
    case 4
        c=4;
    case 5
        c=5;
    case 6
        c=6;
    case 7
        c=7;
    case 8
        c=8;
    case 9
        c=9;
    case 10
        c=10;
end









function mfoc_Callback(hObject, eventdata, handles)
% hObject    handle to mfoctext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc as text
%        str2double(get(hObject,'String')) returns contents of mfoc as a double


% --- Executes during object creation, after setting all properties.
function mfoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function mfoc1_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc1 as text
%        str2double(get(hObject,'String')) returns contents of mfoc1 as a double


% --- Executes during object creation, after setting all properties.
function mfoc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc2_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc2 as text
%        str2double(get(hObject,'String')) returns contents of mfoc2 as a double


% --- Executes during object creation, after setting all properties.
function mfoc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc3_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc3 as text
%        str2double(get(hObject,'String')) returns contents of mfoc3 as a double


% --- Executes during object creation, after setting all properties.
function mfoc3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc4_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc4 as text
%        str2double(get(hObject,'String')) returns contents of mfoc4 as a double


% --- Executes during object creation, after setting all properties.
function mfoc4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc5_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc5 as text
%        str2double(get(hObject,'String')) returns contents of mfoc5 as a double


% --- Executes during object creation, after setting all properties.
function mfoc5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc6_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc6 as text
%        str2double(get(hObject,'String')) returns contents of mfoc6 as a double


% --- Executes during object creation, after setting all properties.
function mfoc6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc7_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc7 as text
%        str2double(get(hObject,'String')) returns contents of mfoc7 as a double


% --- Executes during object creation, after setting all properties.
function mfoc7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc8_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc8 as text
%        str2double(get(hObject,'String')) returns contents of mfoc8 as a double


% --- Executes during object creation, after setting all properties.
function mfoc8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc9_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc9 as text
%        str2double(get(hObject,'String')) returns contents of mfoc9 as a double


% --- Executes during object creation, after setting all properties.
function mfoc9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc10_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc10 as text
%        str2double(get(hObject,'String')) returns contents of mfoc10 as a double


% --- Executes during object creation, after setting all properties.
function mfoc10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in completed.
function completed_Callback(hObject, eventdata, handles)


% mfoc(1)=str2double(get(handles.mfoc1, 'String'));
% mfoc(2)=str2double(get(handles.mfoc2, 'String'));
% mfoc(3)=str2double(get(handles.mfoc3, 'String'));
% mfoc(4)=str2double(get(handles.mfoc4, 'String'));
% mfoc(5)=str2double(get(handles.mfoc5, 'String'));
% mfoc(6)=str2double(get(handles.mfoc6, 'String'));
% mfoc(7)=str2double(get(handles.mfoc7, 'String'));
% mfoc(8)=str2double(get(handles.mfoc8, 'String'));
% mfoc(9)=str2double(get(handles.mfoc9, 'String'));
% mfoc(10)=str2double(get(handles.mfoc10, 'String'));
% 
% noc=get(handles.noc, 'Value');
% for i=1:noc
%     X(i)=mfoc(i);
% end
% 
% X













function not_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to not_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of not_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of not_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function not_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to not_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit128_Callback(hObject, eventdata, handles)
% hObject    handle to edit128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit128 as text
%        str2double(get(hObject,'String')) returns contents of edit128 as a double


% --- Executes during object creation, after setting all properties.
function edit128_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit128 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ftemp_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to ftemp_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftemp_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of ftemp_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function ftemp_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftemp_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in fs_uipanel19.
function fs_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to fs_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fs_uipanel19 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fs_uipanel19


% --- Executes during object creation, after setting all properties.
function fs_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xd_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to xd_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xd_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of xd_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function xd_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xd_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xw_uipanel19_Callback(hObject, eventdata, handles)
% hObject    handle to xw_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xw_uipanel19 as text
%        str2double(get(hObject,'String')) returns contents of xw_uipanel19 as a double


% --- Executes during object creation, after setting all properties.
function xw_uipanel19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xw_uipanel19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qr_Callback(hObject, eventdata, handles)
% hObject    handle to qr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qr as text
%        str2double(get(hObject,'String')) returns contents of qr as a double


% --- Executes during object creation, after setting all properties.
function qr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qc_Callback(hObject, eventdata, handles)
% hObject    handle to qc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qc as text
%        str2double(get(hObject,'String')) returns contents of qc as a double


% --- Executes during object creation, after setting all properties.
function qc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit133_Callback(hObject, eventdata, handles)
% hObject    handle to edit133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit133 as text
%        str2double(get(hObject,'String')) returns contents of edit133 as a double


% --- Executes during object creation, after setting all properties.
function edit133_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit133 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit134_Callback(hObject, eventdata, handles)
% hObject    handle to edit134 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit134 as text
%        str2double(get(hObject,'String')) returns contents of edit134 as a double


% --- Executes during object creation, after setting all properties.
function edit134_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit134 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit135_Callback(hObject, eventdata, handles)
% hObject    handle to edit135 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit135 as text
%        str2double(get(hObject,'String')) returns contents of edit135 as a double


% --- Executes during object creation, after setting all properties.
function edit135_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit135 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit136_Callback(hObject, eventdata, handles)
% hObject    handle to edit136 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit136 as text
%        str2double(get(hObject,'String')) returns contents of edit136 as a double


% --- Executes during object creation, after setting all properties.
function edit136_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit136 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit137_Callback(hObject, eventdata, handles)
% hObject    handle to edit137 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit137 as text
%        str2double(get(hObject,'String')) returns contents of edit137 as a double


% --- Executes during object creation, after setting all properties.
function edit137_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit137 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit138_Callback(hObject, eventdata, handles)
% hObject    handle to edit138 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit138 as text
%        str2double(get(hObject,'String')) returns contents of edit138 as a double


% --- Executes during object creation, after setting all properties.
function edit138_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit138 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit139_Callback(hObject, eventdata, handles)
% hObject    handle to edit139 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit139 as text
%        str2double(get(hObject,'String')) returns contents of edit139 as a double


% --- Executes during object creation, after setting all properties.
function edit139_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit139 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit140_Callback(hObject, eventdata, handles)
% hObject    handle to edit140 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit140 as text
%        str2double(get(hObject,'String')) returns contents of edit140 as a double


% --- Executes during object creation, after setting all properties.
function edit140_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit140 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in noc_uipanel23.
function noc_uipanel23_Callback(hObject, eventdata, handles)

set(handles.uipanel25, 'Visible', 'on');


popup_sel_index = get(handles.noc_uipanel23, 'Value');
switch popup_sel_index

    case 1
        
        set(handles.mfoc11, 'Visible', 'off');
        set(handles.mfoc22, 'Visible', 'off');
        set(handles.mfoc33, 'Visible', 'off');
        set(handles.mfoc44, 'Visible', 'off');
        set(handles.mfoc55, 'Visible', 'off');
        set(handles.mfoc66, 'Visible', 'off');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
        set(handles.mfoctext11, 'Visible', 'off');
        set(handles.mfoctext22, 'Visible', 'off');
        set(handles.mfoctext33, 'Visible', 'off');
        set(handles.mfoctext44, 'Visible', 'off');
        set(handles.mfoctext55, 'Visible', 'off');
        set(handles.mfoctext66, 'Visible', 'off');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 2
        
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'off');
        set(handles.mfoc44, 'Visible', 'off');
        set(handles.mfoc55, 'Visible', 'off');
        set(handles.mfoc66, 'Visible', 'off');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
        set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'off');
        set(handles.mfoctext44, 'Visible', 'off');
        set(handles.mfoctext55, 'Visible', 'off');
        set(handles.mfoctext66, 'Visible', 'off');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 3
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'off');
        set(handles.mfoc55, 'Visible', 'off');
        set(handles.mfoc66, 'Visible', 'off');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
         set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'off');
        set(handles.mfoctext55, 'Visible', 'off');
        set(handles.mfoctext66, 'Visible', 'off');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 4
         set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'off');
        set(handles.mfoc66, 'Visible', 'off');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
        set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'off');
        set(handles.mfoctext66, 'Visible', 'off');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 5
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'off');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
        set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'off');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 6
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'on');
        set(handles.mfoc77, 'Visible', 'off');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
         set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'on');
        set(handles.mfoctext77, 'Visible', 'off');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 7
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'on');
        set(handles.mfoc77, 'Visible', 'on');
        set(handles.mfoc88, 'Visible', 'off');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
        set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'on');
        set(handles.mfoctext77, 'Visible', 'on');
        set(handles.mfoctext88, 'Visible', 'off');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 8
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'on');
        set(handles.mfoc77, 'Visible', 'on');
        set(handles.mfoc88, 'Visible', 'on');
        set(handles.mfoc99, 'Visible', 'off');
        set(handles.mfoc1010, 'Visible', 'off');
         set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'on');
        set(handles.mfoctext77, 'Visible', 'on');
        set(handles.mfoctext88, 'Visible', 'on');
        set(handles.mfoctext99, 'Visible', 'off');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 9
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'on');
        set(handles.mfoc77, 'Visible', 'on');
        set(handles.mfoc88, 'Visible', 'on');
        set(handles.mfoc99, 'Visible', 'on');
        set(handles.mfoc1010, 'Visible', 'off');
         set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'on');
        set(handles.mfoctext77, 'Visible', 'on');
        set(handles.mfoctext88, 'Visible', 'on');
        set(handles.mfoctext99, 'Visible', 'on');
        set(handles.mfoctext1010, 'Visible', 'off');
    case 10
        set(handles.mfoc11, 'Visible', 'on');
        set(handles.mfoc22, 'Visible', 'on');
        set(handles.mfoc33, 'Visible', 'on');
        set(handles.mfoc44, 'Visible', 'on');
        set(handles.mfoc55, 'Visible', 'on');
        set(handles.mfoc66, 'Visible', 'on');
        set(handles.mfoc77, 'Visible', 'on');
        set(handles.mfoc88, 'Visible', 'on');
        set(handles.mfoc99, 'Visible', 'on');
        set(handles.mfoc1010, 'Visible', 'on');
        set(handles.mfoctext11, 'Visible', 'on');
        set(handles.mfoctext22, 'Visible', 'on');
        set(handles.mfoctext33, 'Visible', 'on');
        set(handles.mfoctext44, 'Visible', 'on');
        set(handles.mfoctext55, 'Visible', 'on');
        set(handles.mfoctext66, 'Visible', 'on');
        set(handles.mfoctext77, 'Visible', 'on');
        set(handles.mfoctext88, 'Visible', 'on');
        set(handles.mfoctext99, 'Visible', 'on');
        set(handles.mfoctext1010, 'Visible', 'on');

end



% --- Executes during object creation, after setting all properties.
function noc_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noc_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit141_Callback(hObject, eventdata, handles)
% hObject    handle to edit141 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit141 as text
%        str2double(get(hObject,'String')) returns contents of edit141 as a double


% --- Executes during object creation, after setting all properties.
function edit141_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit141 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fr_uipanel23_Callback(hObject, eventdata, handles)
% hObject    handle to fr_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_uipanel23 as text
%        str2double(get(hObject,'String')) returns contents of fr_uipanel23 as a double


% --- Executes during object creation, after setting all properties.
function fr_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_uipanel23_Callback(hObject, eventdata, handles)
% hObject    handle to d_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_uipanel23 as text
%        str2double(get(hObject,'String')) returns contents of d_uipanel23 as a double


% --- Executes during object creation, after setting all properties.
function d_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_uipanel23_Callback(hObject, eventdata, handles)
% hObject    handle to r_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_uipanel23 as text
%        str2double(get(hObject,'String')) returns contents of r_uipanel23 as a double


% --- Executes during object creation, after setting all properties.
function r_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit147_Callback(hObject, eventdata, handles)
% hObject    handle to edit147 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit147 as text
%        str2double(get(hObject,'String')) returns contents of edit147 as a double


% --- Executes during object creation, after setting all properties.
function edit147_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit147 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit148_Callback(hObject, eventdata, handles)
% hObject    handle to edit148 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit148 as text
%        str2double(get(hObject,'String')) returns contents of edit148 as a double


% --- Executes during object creation, after setting all properties.
function edit148_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit148 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit149_Callback(hObject, eventdata, handles)
% hObject    handle to edit149 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit149 as text
%        str2double(get(hObject,'String')) returns contents of edit149 as a double


% --- Executes during object creation, after setting all properties.
function edit149_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit149 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit150_Callback(hObject, eventdata, handles)
% hObject    handle to edit150 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit150 as text
%        str2double(get(hObject,'String')) returns contents of edit150 as a double


% --- Executes during object creation, after setting all properties.
function edit150_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit150 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit151_Callback(hObject, eventdata, handles)
% hObject    handle to edit151 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit151 as text
%        str2double(get(hObject,'String')) returns contents of edit151 as a double


% --- Executes during object creation, after setting all properties.
function edit151_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit151 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit152_Callback(hObject, eventdata, handles)
% hObject    handle to edit152 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit152 as text
%        str2double(get(hObject,'String')) returns contents of edit152 as a double


% --- Executes during object creation, after setting all properties.
function edit152_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit152 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit153_Callback(hObject, eventdata, handles)
% hObject    handle to edit153 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit153 as text
%        str2double(get(hObject,'String')) returns contents of edit153 as a double


% --- Executes during object creation, after setting all properties.
function edit153_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit153 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit154_Callback(hObject, eventdata, handles)
% hObject    handle to edit154 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit154 as text
%        str2double(get(hObject,'String')) returns contents of edit154 as a double


% --- Executes during object creation, after setting all properties.
function edit154_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit154 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit155_Callback(hObject, eventdata, handles)
% hObject    handle to edit155 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit155 as text
%        str2double(get(hObject,'String')) returns contents of edit155 as a double


% --- Executes during object creation, after setting all properties.
function edit155_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit155 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit156_Callback(hObject, eventdata, handles)
% hObject    handle to edit156 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit156 as text
%        str2double(get(hObject,'String')) returns contents of edit156 as a double


% --- Executes during object creation, after setting all properties.
function edit156_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit156 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit157_Callback(hObject, eventdata, handles)
% hObject    handle to edit157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit157 as text
%        str2double(get(hObject,'String')) returns contents of edit157 as a double


% --- Executes during object creation, after setting all properties.
function edit157_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit157 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit158_Callback(hObject, eventdata, handles)
% hObject    handle to edit158 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit158 as text
%        str2double(get(hObject,'String')) returns contents of edit158 as a double


% --- Executes during object creation, after setting all properties.
function edit158_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit158 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit159_Callback(hObject, eventdata, handles)
% hObject    handle to edit159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit159 as text
%        str2double(get(hObject,'String')) returns contents of edit159 as a double


% --- Executes during object creation, after setting all properties.
function edit159_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit159 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit160_Callback(hObject, eventdata, handles)
% hObject    handle to edit160 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit160 as text
%        str2double(get(hObject,'String')) returns contents of edit160 as a double


% --- Executes during object creation, after setting all properties.
function edit160_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit160 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit161_Callback(hObject, eventdata, handles)
% hObject    handle to edit161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit161 as text
%        str2double(get(hObject,'String')) returns contents of edit161 as a double


% --- Executes during object creation, after setting all properties.
function edit161_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit161 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit162_Callback(hObject, eventdata, handles)
% hObject    handle to edit162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit162 as text
%        str2double(get(hObject,'String')) returns contents of edit162 as a double


% --- Executes during object creation, after setting all properties.
function edit162_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit162 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit163_Callback(hObject, eventdata, handles)
% hObject    handle to edit163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit163 as text
%        str2double(get(hObject,'String')) returns contents of edit163 as a double


% --- Executes during object creation, after setting all properties.
function edit163_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit163 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit164_Callback(hObject, eventdata, handles)
% hObject    handle to edit164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit164 as text
%        str2double(get(hObject,'String')) returns contents of edit164 as a double


% --- Executes during object creation, after setting all properties.
function edit164_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit164 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fs_uipanel23.
function fs_uipanel23_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function fs_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in dp_uipanel23.
function dp_uipanel23_Callback(hObject, eventdata, handles)

set(handles.fr_uipanel23 , 'String',  100  );
set(handles.d_uipanel23 , 'String',  50  );
set(handles.r_uipanel23 , 'String',  1  );
set(handles.ftemp_uipanel23 , 'String',  400  );
set(handles.bldlupper , 'String',  5e-4  );
set(handles.bhdhlower , 'String',  3.1  );


set(handles.noc_uipanel23 , 'Value',  3  );

set(handles.uipanel25, 'Visible', 'on');
set(handles.mfoc11 , 'String', 1/3  );
set(handles.mfoc22 , 'String', 1/3  );
set(handles.mfoc33 , 'String', 1/3  );
set(handles.mfoc11, 'Visible', 'on');
set(handles.mfoc22, 'Visible', 'on');
set(handles.mfoc33, 'Visible', 'on');
set(handles.mfoctext11, 'Visible', 'on');
set(handles.mfoctext22, 'Visible', 'on');
set(handles.mfoctext33, 'Visible', 'on');

set(handles.mfoc44, 'Visible', 'off');
set(handles.mfoc55, 'Visible', 'off');
set(handles.mfoc66, 'Visible', 'off');
set(handles.mfoc77, 'Visible', 'off');
set(handles.mfoc88, 'Visible', 'off');
set(handles.mfoc99, 'Visible', 'off');
set(handles.mfoc1010, 'Visible', 'off');
set(handles.mfoctext44, 'Visible', 'off');
set(handles.mfoctext55, 'Visible', 'off');
set(handles.mfoctext66, 'Visible', 'off');
set(handles.mfoctext77, 'Visible', 'off');
set(handles.mfoctext88, 'Visible', 'off');
set(handles.mfoctext99, 'Visible', 'off');
set(handles.mfoctext1010, 'Visible', 'off');

set(handles.fs_uipanel23 , 'Value', 4  );

set(handles.pw, 'Visible', 'off');

% --- Executes on button press in back_uipanel23.
function pushbutton59_Callback(hObject, eventdata, handles)
% hObject    handle to back_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit165_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc1 as text
%        str2double(get(hObject,'String')) returns contents of mfoc1 as a double


% --- Executes during object creation, after setting all properties.
function edit165_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit166_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc2 as text
%        str2double(get(hObject,'String')) returns contents of mfoc2 as a double


% --- Executes during object creation, after setting all properties.
function edit166_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit167_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc3 as text
%        str2double(get(hObject,'String')) returns contents of mfoc3 as a double


% --- Executes during object creation, after setting all properties.
function edit167_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit168_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc4 as text
%        str2double(get(hObject,'String')) returns contents of mfoc4 as a double


% --- Executes during object creation, after setting all properties.
function edit168_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit169_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc5 as text
%        str2double(get(hObject,'String')) returns contents of mfoc5 as a double


% --- Executes during object creation, after setting all properties.
function edit169_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit170_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc6 as text
%        str2double(get(hObject,'String')) returns contents of mfoc6 as a double


% --- Executes during object creation, after setting all properties.
function edit170_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit171_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc7 as text
%        str2double(get(hObject,'String')) returns contents of mfoc7 as a double


% --- Executes during object creation, after setting all properties.
function edit171_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit172_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc8 as text
%        str2double(get(hObject,'String')) returns contents of mfoc8 as a double


% --- Executes during object creation, after setting all properties.
function edit172_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit173_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc9 as text
%        str2double(get(hObject,'String')) returns contents of mfoc9 as a double


% --- Executes during object creation, after setting all properties.
function edit173_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit174_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc10 as text
%        str2double(get(hObject,'String')) returns contents of mfoc10 as a double


% --- Executes during object creation, after setting all properties.
function edit174_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in back_uipanel7.
function back_uipanel7_Callback(hObject, eventdata, handles)


set(handles.figure1, 'Position', [95,25,65,15.5]);
set(handles.uipanel7, 'Visible', 'off');
set(handles.uipanel5, 'Visible', 'on');
set(handles.uipanel5, 'Position', [2.4,0.692,60.6,14.231]);
set(handles.uipanel9, 'Visible', 'off');
set(handles.uipanel10, 'Visible', 'off');
set(handles.uipanel11, 'Visible', 'off');
set(handles.figure1, 'Name','Process Optimizer');
set(handles.uipanel13, 'Visible', 'off');
set(handles.uipanel16, 'Visible', 'off');
set(handles.uipanel17, 'Visible', 'off');
set(handles.back_uipanel7, 'Visible', 'off');
set(handles.uipanel26, 'Visible', 'off');






% --- Executes on button press in pushbutton62.
function pushbutton62_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dp_uipanel26.
function dp_uipanel26_Callback(hObject, eventdata, handles)


ra=get(handles.ra,'Value');
ga=get(handles.ga,'Value');
pso=get(handles.pso,'Value');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ra==1
    % Rain (defalt parameters)
    set(handles.nov_uipanel9,'String','3');
    set(handles.pp_uipanel9,'String','50');
    set(handles.nop_uipanel9,'String','11');
    set(handles.noi_uipanel9,'String','1');
    set(handles.norp_uipanel9,'String','11');
    set(handles.nord_uipanel9,'String','9');
    set(handles.mint_uipanel9,'String','4');
    set(handles.maxt_uipanel9,'String','100');

    set(handles.uipanel10, 'Visible', 'off');
    set(handles.uipanel11, 'Visible', 'off');

    set(handles.uipanel9, 'Visible', 'on');
    set(handles.uipanel9, 'Position', [33,1.2,39.6,25.846]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ga==1
    % GA (defalt parameters)
    set(handles.nov_uipanel10,'String','3');
    set(handles.nop_uipanel10,'String','100');
    set(handles.nog_uipanel10,'String','20');
    set(handles.nomcg_uipanel10,'String','20');
    set(handles.nomcr_uipanel10,'String','20');
    set(handles.noec_uipanel10,'String','2');
    set(handles.mint_uipanel10,'String','4');
    set(handles.maxt_uipanel10,'String','100');

    set(handles.uipanel9, 'Visible', 'off');
    set(handles.uipanel11, 'Visible', 'off');

    set(handles.uipanel10, 'Visible', 'on');
    set(handles.uipanel10, 'Position', [33,1.2,39.6,25.846]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pso==1
    % PSO (defalt parameters)
    set(handles.nov_uipanel11,'String','3');
    set(handles.pa_uipanel11,'String','set');
    set(handles.pb_uipanel11,'String','set');
    set(handles.pc_uipanel11,'String','1');
    set(handles.pd_uipanel11,'String','1');
    set(handles.ivov_uipanel11,'String','-0.1');
    set(handles.ivox_uipanel11,'String','1');
    set(handles.noi_uipanel11,'String','500');
    set(handles.mint_uipanel11,'String','4');
    set(handles.maxt_uipanel11,'String','100');


    set(handles.uipanel9, 'Visible', 'off');
    set(handles.uipanel10, 'Visible', 'off');
    %
    set(handles.uipanel11, 'Visible', 'on');
    set(handles.uipanel11, 'Position', [33,1.2,39.6,25.846]);
end







function mfoc11_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc11 as text
%        str2double(get(hObject,'String')) returns contents of mfoc11 as a double


% --- Executes during object creation, after setting all properties.
function mfoc11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc22_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc22 as text
%        str2double(get(hObject,'String')) returns contents of mfoc22 as a double


% --- Executes during object creation, after setting all properties.
function mfoc22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc33_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc33 as text
%        str2double(get(hObject,'String')) returns contents of mfoc33 as a double


% --- Executes during object creation, after setting all properties.
function mfoc33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc44_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc44 as text
%        str2double(get(hObject,'String')) returns contents of mfoc44 as a double


% --- Executes during object creation, after setting all properties.
function mfoc44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc55_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc55 as text
%        str2double(get(hObject,'String')) returns contents of mfoc55 as a double


% --- Executes during object creation, after setting all properties.
function mfoc55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc66_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc66 as text
%        str2double(get(hObject,'String')) returns contents of mfoc66 as a double


% --- Executes during object creation, after setting all properties.
function mfoc66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc77_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc77 as text
%        str2double(get(hObject,'String')) returns contents of mfoc77 as a double


% --- Executes during object creation, after setting all properties.
function mfoc77_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc88_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc88 as text
%        str2double(get(hObject,'String')) returns contents of mfoc88 as a double


% --- Executes during object creation, after setting all properties.
function mfoc88_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc99_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc99 as text
%        str2double(get(hObject,'String')) returns contents of mfoc99 as a double


% --- Executes during object creation, after setting all properties.
function mfoc99_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mfoc1010_Callback(hObject, eventdata, handles)
% hObject    handle to mfoc1010 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfoc1010 as text
%        str2double(get(hObject,'String')) returns contents of mfoc1010 as a double


% --- Executes during object creation, after setting all properties.
function mfoc1010_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfoc1010 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit185_Callback(hObject, eventdata, handles)
% hObject    handle to edit185 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit185 as text
%        str2double(get(hObject,'String')) returns contents of edit185 as a double


% --- Executes during object creation, after setting all properties.
function edit185_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit185 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ftemp_uipanel23_Callback(hObject, eventdata, handles)
% hObject    handle to ftemp_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftemp_uipanel23 as text
%        str2double(get(hObject,'String')) returns contents of ftemp_uipanel23 as a double


% --- Executes during object creation, after setting all properties.
function ftemp_uipanel23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftemp_uipanel23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ra.
function ra_Callback(hObject, eventdata, handles)
% hObject    handle to ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ra


% --- Executes on button press in ga.
function ga_Callback(hObject, eventdata, handles)
% hObject    handle to ga (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ga


% --- Executes on button press in pso.
function pso_Callback(hObject, eventdata, handles)
% hObject    handle to pso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pso





function bldlupper_Callback(hObject, eventdata, handles)
% hObject    handle to bldlupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bldlupper as text
%        str2double(get(hObject,'String')) returns contents of bldlupper as a double


% --- Executes during object creation, after setting all properties.
function bldlupper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bldlupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bhdhlower_Callback(hObject, eventdata, handles)
% hObject    handle to bhdhlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bhdhlower as text
%        str2double(get(hObject,'String')) returns contents of bhdhlower as a double


% --- Executes during object creation, after setting all properties.
function bhdhlower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bhdhlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


