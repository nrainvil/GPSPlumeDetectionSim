clear all;

Er = 6371; %km
GPSr = 26600; %km 
pWidth = .4; %km 
rcvr2Vent = 2; %km

el = 30;

theta0 = 90+el;
theta1 = asind(Er*sind(theta0)/GPSr);
theta2 = 180 - theta1 - theta0;
fprintf('A:%0.1f B:%0.1f C:%0.1f\n',theta0,theta1,theta2);

rcvr2GPS = sind(theta2)*GPSr/sind(theta0); %km
thetaP = atand(pWidth/(rcvr2GPS-rcvr2Vent));
tSpace = rcvr2GPS*pWidth/(rcvr2GPS-rcvr2Vent);

fprintf('r2G: %0.0f km thetaP: %0.1f\n',rcvr2GPS,thetaP);
fprintf('pWidth: %0.2f m tSpace: %0.2f m\n',1000*pWidth,1000*tSpace);

h = rcvr2Vent*tand(el);
fprintf('h: %0.2f m\n',1000*h);