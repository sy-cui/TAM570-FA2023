function [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Q,Bl,Rx,Dh);

%%   Evaluate curl-curl term (to be extrapolated)

%    Omega = Lxi*(Dh*V) - Lyi*(U*Dh');
%    curlcurlX =  Bxy.*(Lyi*(Omega*Dh'));
%    curlcurlY = -Bxy.*(Lxi*(Dh*Omega));
     
     rx = Rx(:,:,:,1,1);  % dr/dx
     ry = Rx(:,:,:,1,2);  % dr/dy
     sx = Rx(:,:,:,2,1);  % ds/dx
     sy = Rx(:,:,:,2,2);  % ds/dy 

     Omega = tensor3(1,1,Dh,V).*rx + tensor3(Dh,1,1,V).*sx ...
          - tensor3(1,1,Dh,U).*ry - tensor3(Dh,1,1,U).*sy;
     curlcurlX = ry.*tensor3(1,1,Dh,Omega) + sy.*tensor3(Dh,1,1,Omega);
     curlcurlY = -rx.*tensor3(1,1,Dh,Omega) - sx.*tensor3(Dh,1,1,Omega);
     curlcurlX = Bl.*curlcurlX;
     curlcurlY = Bl.*curlcurlY;

     % Omega = 0*U;
     % curlcurlX = 0*U;
     % curlcurlY = 0*U;
