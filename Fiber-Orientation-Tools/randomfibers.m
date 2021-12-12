function p = randomfibers(N)
% function p = randomfibers(N);
% Generate N 3-D randomly oriented fibers
% p is a (3xn) matrix of unit orientation vector components

p = zeros(3,N);  % Storage for the results

% Generate fibers, by angle and by p vector
phi = 2*pi* rand(N, 1);         % angles phi of the fibers
z  = rand(N,1);                 % Ahmed's z
theta = acos(1-2*z);            % Chuck's theta = q(z)
p(1,:) = sin(theta).*cos(phi);  % p1 vector components
p(2,:) = sin(theta).*sin(phi);  % p2 vector components
p(3,:) = cos(theta);            % p3 vector components (= Ahmed's z)

% Plot the family of fibers thus generated (N = 500 gives reasonable results)
% plot3(p(1,:),p(2,:),p(3,:),'r*'); axis equal

return;