datacorr = csvread('2Dperiod_SpacialCorr.csv');
hold off;

%% 横軸リスケーリング
hold off;
hold on;
datax = 1:149;
for i = 1:149
    datax(1,i) = i;
end
x = zeros(5,2);
% x = zeros(5,1);
t = 0;
for i = [6,11,31,51,101]
    t = t+1;
    name = strcat('t=',num2str(1000*(i-1)));
    plot(datax/sqrt(1000*(i-1)),datacorr(i,:)*log(1000*(i-1)/300), 'DisplayName', name);
    ydata = datacorr(i,:)*log(1000*(i-1)/300);
%     func = @(x,datax)expint(x(1)*datax.^2/(2*1000*(i-1)));
    func = @(x,datax)x(2)*expint(x(1)*datax.^2/(2*1000*(i-1)));
    x0 = [1 1];
    x(t,:) = lsqcurvefit(func,x0,datax,ydata);
end
legend('show');


%% 相関関数プロット
datax = 1:149;
for i = 1:149
    datax(1,i) = i;
end
axes('Position',[.7 .7 .2 .2]);
box on;
for i = [6,11,31,51,101]
    name = strcat('t=',num2str(1000*(i-1)));
    plot(datax,datacorr(i,:), 'DisplayName', name);
    hold on
end
legend('show');


%% 指数積分関数のフィッティング
datax = 1:500;
cy = 1:500;
t1 = 100000;
cons = 3;
for i = 1:500
    cy(1,i) = x(5,2)*expint(x(5,1)*cy(1,i)^2/(2*t1));
end
plot(datax/sqrt(t1),cy,'DisplayName','E.i. function');
