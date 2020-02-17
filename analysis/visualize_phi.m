mov = VideoWriter('movc.avi');
mov.FrameRate = 60;
open(mov);

amp = zeros(1,5000);
for t = 1:5000
    name = strcat("phi",num2str(t*2),".csv");
    data = csvread(name);
    
    data(:,151) = [];
    surf(data);
    zlim([0 1]);
    title(['t=',num2str(2*t)],'Fontsize',20);
    F(t) = getframe(gcf);
    writeVideo(mov,F(t));
    
%     for x = 1:150
%         for y = 1:150
%             amp(1,t) = amp(1,t)+abs((data(x,y)-0.5)*2);
%         end
%     end
%     amp(1,t) = amp(1,t)/(150*150);
end

close(mov);

% ampx = 1:5001;
% plot(ampx,amp);