% %Brownian motion
% clf;
% n=20;
% s=0.02;
% x = rand(n,1)-0.5;
% y = rand(n,1)-0.5;
% h = plot(x,y,'.');
% axis([-1 1 -1 1])
% axis square
% grid off
% set(h,'EraseMode','xor','MarkerSize',18)
% grid on;
% title('Press Ctl-C to stop');
% while 1
% drawnow
% x = x + s*randn(n,1);
% y = y + s*randn(n,1);
% set(h,'XData',x,'YData',y)
% end

function o=brownian(n)
% beta=1.5;
% sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% u=randn(1,d)*sigma;
% v=randn(1,d);
% step=u./abs(v).^(1/beta);
% o=step;
s=0.02;
o=s*randn(n,1);
end