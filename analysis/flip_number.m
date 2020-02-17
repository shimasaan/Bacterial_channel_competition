datand = csvread('2Dperiod_n_distribution.csv');
datan = csvread('2Dperiod_n.csv');

%%
datax = 1:3000;
hold on;
for i = [2,6,11,31,51,101]
    name = strcat('t=',num2str(1000*(i-1)));
    plot(datax, datand(i,:), 'DisplayName', name);
end
legend('show');

%%
clx = 1:100001;
for i = 1:100001
    datan(i,1) = datan(i,1)*log(i);
end
axes('Position',[.7 .7 .2 .2]);
box on;
plot(clx,datan);