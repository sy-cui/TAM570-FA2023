function [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Q,Bl,Rx,Dh);

%%   Evaluate curl-curl term (to be extrapolated)

%    Omega = Lxi*(Dh*V) - Lyi*(U*Dh');
%    curlcurlX =  Bxy.*(Lyi*(Omega*Dh'));
%    curlcurlY = -Bxy.*(Lxi*(Dh*Omega));
     
     rx = Rx(:,:,:,1,1);  % dr/dx
     ry = Rx(:,:,:,1,2);  % dr/dy
     sx = Rx(:,:,:,2,1);  % ds/dx
     sy = Rx(:,:,:,2,2);  % ds/dy 

     ur = V.*rx - U.*ry;
     us = V.*sx - U.*sy;

     Omega = tensor3(1,1,Dh,ur) + tensor3(Dh,1,1,us);
     curlcurlX = tensor3(1,1,Dh,ry.*Omega) + tensor3(Dh,1,1,sy.*Omega);
     curlcurlY = -tensor3(1,1,Dh,rx.*Omega) - tensor3(Dh,1,1,sx.*Omega);
     curlcurlX = Bl.*curlcurlX;
     curlcurlY = Bl.*curlcurlY;

     % Omega = 0*U;
     % curlcurlX = 0*U;
     % curlcurlY = 0*U;
