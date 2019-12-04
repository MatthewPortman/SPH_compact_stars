system('make');
system('./sph_test');

% Using parameters from the .F90 code
dt = 0.04;
it_t = 10;
tmax = it_t/dt;
tgrid = 0 : tmax;
nt = length(tgrid);
maxn = 1000;
dx = maxn/2;
xgrid = -1 : dx : 1;

pos = importdata('./1000particles/pos_1000.dat', ' ');
ut = importdata('./1000particles/ut_1000.dat',' ');

% Initializing and reshaping
posi = zeros(maxn,2,tmax);

ut = reshape(ut,maxn,tmax);
posi(1:maxn,:,1) = pos(1:maxn,:);

% Adjusting the position array to be used in my graph.
for i = 1:tmax-1
    posi(1:maxn,:,i+1) = pos(i*maxn+1:(i+1)*maxn,:);
end

% Plotting and giffing
str = strjoin({'Toy_Star_2D_',num2str(maxn),'.gif'});

     for jj = 1:tmax
   
        drawnow
     
        % Plotting the results.
        scatter(posi(:,1,jj), posi(:,2,jj),1,ut(:,jj))
        %hold on;
        colorbar
        title(['Toy Star 2D \newline t = ' num2str(jj*dt) ', dt = ' num2str(0.04)])
        xlabel('x')
        ylabel('y')
        zlabel('ut')
        axis([-1.5 1.5 -1.5 1.5])  
        grid on
        
        if jj == 1
          gif(str,'frame',gcf)
        else
          gif
        end
        
%        hold off;
     end