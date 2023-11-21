function anrm = a_norm(Ul,h0,h1,Bl,Grr,Grs,Gss,Dh);

Wl = axl(Ul,h0,h1,Bl,Grr,Grs,Gss,Dh);

anrm = sum(sum(sum(Wl.*Ul)));
anrm = sqrt(anrm);

