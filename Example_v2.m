% Test Example
%
%  Example 1 computed using angles
%
% We consider the Earth magnetic field measurements made by one of the three satellites of the European Space Agency Swarm
% mission from April 21 to 22, 2004. We thank the European Space Agency (ESA) that supports the Swarm mission.
% Swarm data can be accessed online at http://earth.esa.int/swarm.

clear all
clc
savefig=0;
load dataset_preprocessed
%% Signal
FIG=figure;
plot(time,H,'linewidth',2)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 36);
x1=xlabel('time') ; 
set(x1,'Interpreter','latex')
y1=ylabel('H') ; 
set(y1,'Interpreter','latex')
axis tight
if savefig==1
saveas(FIG,'Signal_H', 'fig')
saveas(FIG,'Signal_H', 'epsc')
saveas(FIG,'Signal_H', 'png')
end
%
FIG2=figure
plot(time,D,'linewidth',2)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'fontsize', 36);
x2=xlabel('time') ; 
set(x2,'Interpreter','latex')
y2=ylabel('D') ; 
set(y2,'Interpreter','latex')
axis tight
if savefig==1
saveas(FIG2,'Signal_D', 'fig')
saveas(FIG2,'Signal_D', 'epsc')
saveas(FIG2,'Signal_D', 'png')
end
%
FIG3=figure
plot3(time,H,D,'linewidth',2)
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x3=xlabel('time') ; 
set(x3,'Interpreter','latex')
y3=ylabel('H') ; 
set(y3,'Interpreter','latex')
z3=zlabel('D') ; 
set(z3,'Interpreter','latex')
set(gca,'fontsize', 36);
view([-5 30])
if savefig==1
saveas(FIG3,'Signal', 'fig')
saveas(FIG3,'Signal', 'epsc')
saveas(FIG3,'Signal', 'png')
end

% %% CASE 2D + time
% S=[H,D];
% opts=Settings_FIF_v3('ExtPoints',5,'verbose',0);
% tic
% [IMF_v7,logM_v7] = MvFIF_v7(S,opts);
% toc

%% CASE 3D + time
S=[H,D,Z];
opts=Settings_IF_v1('IF.ExtPoints',5,'verbose',0);
tic
[IMF_v7,logM_v7] = MvFIF_v7(S,opts);
toc

%% Plotting of the first two dimensions 

FIG4=figure;
for ii=1:size(IMF_v7{1},2)-1
    
    plot3(time,IMF_v7{1}(:,ii),IMF_v7{2}(:,ii),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{IMF}_{' num2str(ii) '}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_' num2str(ii)], 'png')
    end
    pause(0.1)
end

ii=size(IMF_v7{1},2)
plot3(time,IMF_v7{1}(:,ii),IMF_v7{2}(:,ii),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{Trend}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_' num2str(ii)], 'png')
    end
    pause(0.1)

    %% CASE 2D + time with minimal number of IMFs (max mask length selected via alpha paramter set to 100)
S=[H,D];
opts=Settings_FIF_v3('ExtPoints',5,'verbose',0,'alpha',100);
tic
[IMF_v7,logM_v7] = MvFIF_v7(S,opts);
toc

%%

FIG4=figure;
for ii=1:13
    
    plot3(time,IMF_v7{1}(:,ii),IMF_v7{2}(:,ii),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{IMF}_{' num2str(ii) '}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)
end

ii=14;
plot3(time,sum(IMF_v7{1}(:,14:17),2),sum(IMF_v7{2}(:,14:17),2),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{IMF}_{' num2str(ii) '}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)
    
    ii=15
    plot3(time,sum(IMF_v7{1}(:,18),2),sum(IMF_v7{2}(:,18),2),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{Trend}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)

%% CASE 3D + time with minimal number of IMFs (max mask length selected via alpha paramter set to 100)
S=[H,D,Z];
opts=Settings_FIF_v3('ExtPoints',5,'verbose',0,'alpha',100);
tic
[IMF_v7,logM_v7] = MvFIF_v7(S,opts);
toc

%%

FIG4=figure;
for ii=1:10
    
    plot3(time,IMF_v7{1}(:,ii),IMF_v7{2}(:,ii),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{IMF}_{' num2str(ii) '}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)
end

ii=11;
plot3(time,sum(IMF_v7{1}(:,11:15),2),sum(IMF_v7{2}(:,11:15),2),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{IMF}_{' num2str(ii) '}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)
    
    ii=12;
    plot3(time,sum(IMF_v7{1}(:,16),2),sum(IMF_v7{2}(:,16),2),'linewidth',2)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    x3=xlabel('time') ;
    set(x3,'Interpreter','latex')
    y3=ylabel('H') ;
    set(y3,'Interpreter','latex')
    z3=zlabel('D') ;
    set(z3,'Interpreter','latex')
    set(gca,'fontsize', 36);
    view([-5 30])
    th=title(['$\textrm{Trend}$']);
    set(th,'Interpreter','latex')
    if savefig==1
    saveas(FIG4,['IMF_max_' num2str(ii)], 'fig')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'epsc')
    saveas(FIG4,['IMF_max_' num2str(ii)], 'png')
    end
    pause(0.1)

