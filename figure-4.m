% equation 71 from Denton et al.
x1 = 0.1;
x2 = 0.3;
x3 = 0.2;
x4 = 0.4;

p0 = x1+x2;
q0 = x1+x3;

del0 = x1*x4-x2*x3;
beta_1_vec = 0:0.001:0.5;
beta_2_vec = [0.01 0.1 0.2 0.3 0.5];

q_eq = zeros(numel(beta_2_vec),numel(beta_1_vec));

for j = 1:numel(beta_2_vec)
    beta_2 = beta_2_vec(j);
    for i = 1:numel(beta_1_vec)
        beta_1 = beta_1_vec(i);
        q_eq(j,i) = q0+(((beta_2-beta_1)*del0)/(beta_1*(1-p0)+beta_2*p0));
    end
end

figure
hold on
plot(beta_1_vec,q0-q_eq,'LineWidth',1.5)
yline(0,'k--','LineWidth',2)
ylabel('Change in q')
xlabel('Beta_1')
legend('\beta_2=0.01','\beta_2=0.1','\beta_2=0.2','\beta_2=0.3','\beta_2=0.5')
hold on 
scatter(beta_2_vec(1),0,'o','filled','HandleVisibility','off')
scatter(beta_2_vec(2),0,'o','filled','HandleVisibility','off')
scatter(beta_2_vec(3),0,'o','filled','HandleVisibility','off')
scatter(beta_2_vec(4),0,'o','filled','HandleVisibility','off')
scatter(beta_2_vec(5),0,'o','filled','HandleVisibility','off')
box on