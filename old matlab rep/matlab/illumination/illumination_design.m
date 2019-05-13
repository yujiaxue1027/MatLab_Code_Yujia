%% illumination design
for theta = 10:5:45 % in degree
for h = [2  4  6]% in mm
for led_rows = [1 2]
for led_cols = [3 4]
L = 10;
dim = 100;
delta = L/dim;
x = -delta*dim/2:delta:delta*dim/2-delta;
y = -delta*dim/2:delta:delta*dim/2-delta;
A = [2,2,10];
B = [-2,2,10];
C = [-2,-2,10];
D = [2,-2,10];
ext = h/tan(theta*pi/180);
AA = A + [ext,ext,-h];
BB = B + [-ext,ext,-h];
CC = C + [-ext,-ext,-h];
DD = D + [ext,-ext,-h];
distribution = zeros(dim,dim);
for i = 1:dim
    for j = 1:dim
        pos = [x(i),y(j),0];
        zunit = [0 0 1]';
        leds = zeros(3,led_rows,led_cols);
        tl = A;
        tr = B;
        bl = AA;
        br = BB;
        vec1 = tr-tl;
        vec2 = bl-tl;
        normal = cross(vec1,vec2);
        normal = normal./norm(normal);
        for ii = 1:led_rows
            for jj = 1:led_cols
                leds(:,ii,jj) = ((led_rows+1-ii)*(led_cols+1-jj)*tl+...
                    (led_rows+1-ii)*jj*tr+...
                    ii*(led_cols+1-jj)*bl+...
                    ii*jj*br)/((led_rows+1)*(led_cols+1));
            end
        end
        for ii = 1:led_rows
            for jj = 1:led_cols
                led_pos = leds(:,ii,jj);
                vec1 = pos(:) - led_pos(:);
                vec1 = vec1./norm(vec1);
                tmp = (dot(vec1,normal));
                vec2 = led_pos(:) - pos(:);
                dist = norm(vec2);
                vec2 = vec2./norm(vec2);                
                tmp = tmp*(dot(vec2,zunit));
                tmp = tmp/((dist)^2);
                distribution(i,j) = distribution(i,j) + tmp;
            end
         end
    end
end
distribution = distribution+rot90(distribution,1)+rot90(distribution,2)+rot90(distribution,3);
distribution = distribution./led_rows./led_cols;
sigma = std(distribution(:));
figure,surf(distribution),zlim([0,0.05]),title(['theta: ',num2str(theta),', h = ',num2str(h),' mm, led: ',...
    num2str(led_rows),' rows, ',num2str(led_cols),' cols, std: ',num2str(sigma)])
end
end
end
end











