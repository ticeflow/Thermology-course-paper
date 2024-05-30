% 2D ising model 模拟
tic
%% MC模拟
% 设置初始参数
L=40;
T=[linspace(0.4,0.8,4),linspace(0.8,1.2,16),linspace(1.2,1.6,4)]*2.27;%设置温度梯度
J=1;  %耦合强度
MCS=1e4;%Monte Carlo步数
step=MCS*L^2; %总步数
H=0; %外场
%初始化
average_spin=zeros(1,step);
sample_spin=zeros(24,MCS,L,L);
E_1=zeros(1,24);
E_2=zeros(1,24);
M_1=zeros(1,24);
M_2=zeros(1,24);
C=zeros(1,24);
X=zeros(1,24);
for tt=1:24
    tic
    %迭代
    % 自旋初态随机取值
    spin=2*randi([0,1],L,L)-1;
    for ii=1:step
        spin=MC(spin,T(tt),H,L,J);
        %采样
        if mod(ii-1,L^2)==0
            sample_spin(tt,(ii-1)/L^2+1,:,:)=spin;
        end
    end
    toc
end
%% 计算能量，磁化强度，比热，磁化率
tic
for tt=1:24
    for ii=500:MCS
        E_1(tt)=E_1(tt)+energy(sample_spin(tt,ii,:,:),L,J)/(MCS-500);
        E_2(tt)=E_2(tt)+energy(sample_spin(tt,ii,:,:),L,J)^2/(MCS-500);
        M_1(tt)=M_1(tt)+abs(sum(sample_spin(tt,ii,:,:),'all'))/(MCS-500);
        M_2(tt)=M_2(tt)+abs(sum(sample_spin(tt,ii,:,:),'all'))^2/(MCS-500);
    end
    C(tt)=(E_2(tt)-E_1(tt)^2)/T(tt)^2/L^2;
    X(tt)=(M_2(tt)-M_1(tt)^2)/T(tt)/L^2;
end  
toc

%% 可视化能量
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
plot(T,E_1/2,'-^r','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel("$T$")
ylabel("$E$")

%% 可视化磁化强度
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
plot(T,M_1,'-^b','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel("$T$")
ylabel("$M$")
%% 可视化比热
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
plot(T,C,'-^y','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel("$T$")
ylabel("$C_v$")
%% 可视化磁化率
figure('Color',[1 1 1])
set(0, 'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
plot(T,X,'-^g','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel("$T$")
ylabel("$\chi$")
%%
function spin=MC(spin,T,H,L,J)
    Exp=[exp(-(-8*J)/T),exp(-(-4*J)/T),1,exp(-(4*J)/T),exp(-(8*J)/T)];
    i=randi([1,L]);
    j=randi([1,L]);
    %计算随机翻转自旋产生的能量变化量
    delta_E=2*J*spin(i,j)*(spin(mod(i,L)+1,j)+spin(mod(i-2,L)+1,j)+spin(i,mod(j,L)+1)+spin(i,mod(j-2,L)+1))+2*H*spin(i,j);
    %此处使用周期性边界条件
    if rand()<Exp(delta_E/(4*J)+3)
        spin(i,j)=-spin(i,j);
    end
end

function E=energy(spin,L,J)
    %计算当前能量
    E=0;
    for ii=1:L
        for jj=1:L
            E=E-J*spin(1,1,ii,jj)*(spin(1,1,mod(ii,L)+1,jj)+spin(1,1,mod(ii-2,L)+1,jj)+spin(1,1,ii,mod(jj,L)+1)+spin(1,1,ii,mod(jj-2,L)+1));
        end
    end
end