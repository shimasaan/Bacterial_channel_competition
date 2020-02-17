space = zeros(900,900);
maxarea = 1000;
F = zeros(5,maxarea);
k = 0;
hold off

%%
datal = 1:maxarea;
hold on;
for i = [2,3,7,11,20]
    k = k+1;
    for ii = 0:9
        name = strcat(num2str(ii),'_',num2str(i),'.csv');
        roi = csvread(name);
        roi(:,301) = [];
        space(1:300,1:300) = roi;
        space(1:300,301:600) = roi;
        space(1:300,601:900) = roi;
        space(301:600,1:300) = roi;
        space(301:600,301:600) = roi;
        space(301:600,601:900) = roi;
        space(601:900,1:300) = roi;
        space(601:900,301:600) = roi;
        space(601:900,601:900) = roi;
        im = imbinarize(space);
        S = regionprops(im,'Area');
        S = struct2array(S);
        si = size(S);
        for j = 1:si(1,2)
            area = 1;
            while S(1,j) >= area & area<=maxarea
                F(k,area) = F(k,area)+1;
                area = area+1;
            end
        end
        im = ~im;
        S = regionprops(im,'Area');
        S = struct2array(S);
        si = size(S);
        for j = 1:si(1,2)
            area = 1;
            while S(1,j) >= area & area<=maxarea
                F(k,area) = F(k,area)+1;
                area = area+1;
            end
        end
    end
    ini = F(k,1);
    for j = 1:maxarea
        F(k,j) = F(k,j)/ini;
    end
end

p1=plot(datal,F(1,:),'DisplayName',strcat('t=',num2str(5000)));
p2=plot(datal,F(2,:),'DisplayName',strcat('t=',num2str(10000)));
p3=plot(datal,F(3,:),'DisplayName',strcat('t=',num2str(30000)));
p4=plot(datal,F(4,:),'DisplayName',strcat('t=',num2str(50000)));
p5=plot(datal,F(5,:),'DisplayName',strcat('t=',num2str(100000)));
x = 1:10:1000;
y = x.^(-0.85);
plot(x,y,'DisplayName','$$a=0.85$$');

legend([p1 p2 p3 p4 p5]);

%%
% clx = 1:100001;
% for i = 1:100001
%     datan(i,1) = datan(i,1)*log(i);
% end
% axes('Position',[.7 .7 .2 .2]);
% box on;
% plot(clx,datan);