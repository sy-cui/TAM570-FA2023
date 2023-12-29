function [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Q,Bl,Rx,Dh);

%%   Evaluate curl-curl term (to be extrapolated)

%%   Be sure that the result is multiplied by Bl() before returning.

%    Omega = Lxi*(Dh*V) - Lyi*(U*Dh');     %%% THIS IS FROM OUR monodomain code
%    curlcurlX =  Bxy.*(Lyi*(Omega*Dh'));  %%% THIS IS FROM OUR monodomain code
%    curlcurlY = -Bxy.*(Lxi*(Dh*Omega));   %%% THIS IS FROM OUR monodomain code

     Omega     = 0*U;
     curlcurlX = Bl.*(0*U);
     curlcurlY = Bl.*(0*U);
