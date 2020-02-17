L = 200;

v = zeros(4*3*L*L,3);
f = zeros(3*L*L,4);
c = zeros(3*L*L,1);
ind = 1;

for i = 1:L %top surface
    for j = 1:L
        v(ind,1) = i-1;
        v(ind,2) = j-1;
        v(ind,3) = L;
        v(ind+1,1) = i;
        v(ind+1,2) = j-1;
        v(ind+1,3) = L;
        v(ind+2,1) = i-1;
        v(ind+2,2) = j;
        v(ind+2,3) = L;
        v(ind+3,1) = i;
        v(ind+3,2) = j;
        v(ind+3,3) = L;
        ind = ind+4;
    end
end

for i = 1:L %left surface
    for j = 1:L
        v(ind,1) = j-1;
        v(ind,2) = 0;
        v(ind,3) = L+1-i;
        v(ind+1,1) = j;
        v(ind+1,2) = 0;
        v(ind+1,3) = L+1-i;
        v(ind+2,1) = j-1;
        v(ind+2,2) = 0;
        v(ind+2,3) = L+1-i-1;
        v(ind+3,1) = j;
        v(ind+3,2) = 0;
        v(ind+3,3) = L+1-i-1;
        ind = ind+4;
    end
end

for i = 1:L %right surface
    for j = 1:L
        v(ind,1) = L;
        v(ind,2) = j-1;
        v(ind,3) = L+1-i;
        v(ind+1,1) = L;
        v(ind+1,2) = j;
        v(ind+1,3) = L+1-i;
        v(ind+2,1) = L;
        v(ind+2,2) = j-1;
        v(ind+2,3) = L+1-i-1;
        v(ind+3,1) = L;
        v(ind+3,2) = j;
        v(ind+3,3) = L+1-i-1;
        ind = ind+4;
    end
end

for i = 1:3*L*L
    f(i,1) = (i-1)*4+1;
    f(i,2) = (i-1)*4+2;
    f(i,3) = (i-1)*4+4;
    f(i,4) = (i-1)*4+3;
end

mov = VideoWriter('movc.avi');
mov.FrameRate = 50;
open(mov);

for t = 1:5001
    namet = strcat('t',num2str(t),'.csv');
    namel = strcat('l',num2str(t),'.csv');
    namer = strcat('r',num2str(t),'.csv');
    top = csvread(namet);
    left = csvread(namel);
    right = csvread(namer);

    ind = 1;
    for i = 1:L
        for j = 1:L
            if top(i,j) == 0
                c(ind,1) = 'y';
            elseif top(i,j) == 1
                c(ind,1) = 'm';
            else
                c(ind,1) = 'r';
            end
%           c(ind,1) = top(i,j);
            ind = ind+1;
        end
    end
    for i = 1:L
        for j = 1:L
            if left(i,j) == 0
                c(ind,1) = 'y';
            elseif left(i,j) == 1
                c(ind,1) = 'm';
            else
                c(ind,1) = 'r';
            end
%           c(ind,1) = left(i,j);
            ind = ind+1;
        end
    end
    for i = 1:L
        for j = 1:L
            if right(i,j) == 0
                c(ind,1) = 'y';
            elseif right(i,j) == 1
                c(ind,1) = 'm';
            else
                c(ind,1) = 'r';
            end
%           c(ind,1) = right(i,j);
            ind = ind+1;
        end
    end

    clf;
    patch('Vertices',[0 0 0; L 0 0; L 0 L; 0 0 L],'faces',[1 2 3 4],'EdgeColor','black','FaceColor','none','LineWidth',1.5);
    patch('Vertices',[L 0 0; L 0 L; L L L; L L 0],'faces',[1 2 3 4],'EdgeColor','black','FaceColor','none','LineWidth',1.5);
    patch('Vertices',[L 0 L; L L L; 0 L L; 0 0 L],'faces',[1 2 3 4],'EdgeColor','black','FaceColor','none','LineWidth',1.5);
    patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceVertexCData',c,'Facecolor','flat');

    title(['t=',num2str(t-1)],'Fontsize',20);
    view(45,25);
    F(t) = getframe(gcf);
    writeVideo(mov,F(t));
end

close(mov);
    