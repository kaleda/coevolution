% Denton et al Paper - 4 or 2 panel eigenvalue figures

%% viability selection
% invasion of h 
% clear all

figure

x1 = 0.5;
x2 = 0.5;
x3 = 0;
x4 = 0;

s_vec = [0:0.01:1];
m_vec = [0:0.01:1];
c_vec = 0:0.1:0.5;

c0 = 0;
c3 = 0.9;

d0 = 0;
d1 = 0.2;
d2 = 0.2;
d3 = 1;

for j = 1:numel(c_vec)
    c1 = c_vec(j);
    c2 = c_vec(j);
    
    alpha = d3*x1+d1*x3;
    beta = c3*x1+c2*x2;
    gamma = d2*x2+d0*x4;
    delta = c1*x3+c0*x4;
    eta = d3*x2+d1*x4;
    phi = c1*x1+c0*x2;
    mu = d2*x1+d0*x3;
    nu = c3*x3+c2*x4;
    
    Z = (alpha*beta)+(gamma*delta)+(eta*phi)+(mu*nu);
    k = d1*c1+d2*c2;
    
    LAM = zeros(1,numel(s_vec));
    
    for i = 1:numel(s_vec)
        s = s_vec(i);
        M = (1+s)*(1+s*Z);
        %J = [(1+s)*(c3*x1+k*x2) (1+s)*((c1+c2)*x1-k*x1); 1-c3*x1-k*x2
        %1-(c1+c2)*x1+k*x1];% only with d1+d2=1!!
        J = [(1+s)*((2-d1-d2)*(c3*x1)+x2*(c1+c2-c1*d2-c2*d1)) (1+s)*((c1+c2)*x1-(d1*c1+c2*d2)*x1); (2-d1-d2)-c3*x1-(c2*(1-d1)+c1*(1-d2))*x2 2-d1-d2-(c1*(1-d1)+c2*(1-d2))*x1];
        LAM(i) = max(eig(J./M));
    end
    
    subplot(1,2,1)
    hold on 
    plot(s_vec,LAM,'LineWidth',1.5)
    ylim([0.6,1.8])
    
    box on 
    
end

yline(1,'k--','LineWidth',1.5)
ylabel('Leading eigenvalue, \lambda_1')
xlabel('Viability selection, s')

%% Fertility selection
% invasion of h 

% x1 = 0.5;
% x2 = 0.5;

% c0 = 0;
% c3 = 0.9;
% 
% d0 = 0;
% d1 = 0.5;
% d2 = 0.5;
% d3 = 1;

for j = 1:numel(c_vec)
    c1 = c_vec(j);
    c2 = c_vec(j);
    
    alpha = d3*x1+d1*x3;
    beta = c3*x1+c2*x2;
    gamma = d2*x2+d0*x4;
    delta = c1*x3+c0*x4;
    eta = d3*x2+d1*x4;
    phi = c1*x1+c0*x2;
    mu = d2*x1+d0*x3;
    nu = c3*x3+c2*x4;
    
    Z = (alpha*beta)+(gamma*delta)+(eta*phi)+(mu*nu);
    k = d1*c1+d2*c2;
    
    LAM = zeros(1,numel(s_vec));
    
    
    for i = 1:numel(m_vec)
        m = m_vec(i);
        M = (1 + m)*(1 + m*x1^2);
        
        de3de3 = -c3*(d1+d2)*x1-(c1*d2+c2*d1)*x2+2*c3*x1+(c1+c2)*x2+m*c3*x1*(2-d1-d2); 
        de3de4 = (c1*(1-d2)+c2*(1-d1))*x1; 
        de4de3 = c3*(d1+d2)*x1+(c1*d2+c2*d1)*x2-2*c3*x1-(c1+c2)*x2+(2-d1-d2)+m*(1-c3)*x1*(2-d1-d2); 
        de4de4 = -(c1*(1-d2)+c2*(1-d1))*x1+(2-d1-d2); 
        
        J = [de3de3 de3de4; de4de3 de4de4];
        LAM(i) = max(eig(J./M));
    end
    
    subplot(1,2,2)
    hold on 
    plot(m_vec,LAM,'LineWidth',1.5)
    ylim([0.6,1.8])
    box on 
    
end

yline(1,'k--','LineWidth',1.5)
xlabel('Fertility selection, m')

legend('c = 0', 'c = 0.1', 'c = 0.2', 'c = 0.3', 'c = 0.4', 'c = 0.5')
