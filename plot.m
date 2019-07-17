inc = 0.01;
theta1 = [0:inc/2:pi];
phi1 = [0:inc:2*pi];
[theta,phi]=meshgrid(theta1,phi1);
l = 6;
leg=legendre(l,cos(theta));
[whatever,row,col]=size(leg);
sph_har=reshape(0*leg(1,:,:),row,col);
Y = zeros(2*l+1,row,col);
s = size(leg);
sp_coeff =[
1.01308e-15,;
-1.35308e-16,;
9.71445e-16,;
-7.63278e-17,;
-6.93889e-16,;
1.94289e-16,;
0.207289,;
-5.96745e-16,;
0.671693,;
2.51882e-15,;
-0.548435,;
-5.55112e-16,;
-0.452856,;
]
;
for m = 0: (s(1)-1)
if m==0
leg0=reshape(leg(1,:,:),row,col);
sph_har=sph_har+sp_coeff(m+1+l)*sqrt((2*l+1)/(4*pi))*leg0;
else
legm=reshape(leg(m+1,:,:),row,col);
sph_har=sph_har+sp_coeff(m+1+l)*(-1)^(m)*sqrt(2)*sqrt((2*l+1)/(4*pi))*sqrt(factorial(l-m)/factorial(l+m))*legm.*cos(m*phi);
sph_har=sph_har+sp_coeff(-m+1+l)*(-1)^(m)*sqrt(2)*sqrt((2*l+1)/(4*pi))*sqrt(factorial(l-m)/factorial(l+m))*legm.*sin(m*phi);
end
end
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta); 
surf(x,y,z,(sph_har)); 
shading interp 
axis equal
