function X=bp(Y,Compressed_Mat)

n=length(Y);                          
s=floor(n/4);                         

X_initial=Compressed_Mat'*Y;                

X=l1eq_pd(s,X_initial, Compressed_Mat, [],Y, 1e-3);



function xp = l1eq_pd(cnt,x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)

% Solve
% min_x ||x||_1  s.t.  Ax = b

largescale = isa(A,'function_handle');

if (nargin < 6), pdtol = 1;  end
if (nargin < 7), pdmaxiter = 5;  end  
if (nargin < 8), cgtol = 1;  end
if (nargin < 9), cgmaxiter = 200;  end

N = length(x0);

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(N,1); ones(N,1)];

if (largescale)
  if (norm(A(x0)-b)/norm(b) > cgtol)

    AAt = @(z) A(At(z));
    [w, cgres, cgiter] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)

      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b)/norm(b) > cgtol)
   
    opts.POSDEF = true; opts.SYM = true;
    [w, hcond] = linsolve(A*A', b, opts);
    if (hcond < 1e-14)
     
      xp = x0;
      return;
    end
    x0 = A'*w;
  end  
end
x = x0;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));

fu1 = x - u;
fu2 = -x - u;
lamu1 = -1./fu1;
lamu2 = -1./fu2;
if (largescale)
  v = -A(lamu1-lamu2);
  Atv = At(v);
  rpri = A(x) - b;
else
  v = -A*(lamu1-lamu2);
  Atv = A'*v;
  rpri = A*x - b;
end

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
done = (sdg < pdtol) | (pditer >= pdmaxiter);
times=0;
while ((~done) && (times<cnt))  
                                
                                
  pditer = pditer + 1;
  times=times+1;
  
  w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
  w3 = -rpri;
  
  sig1 = -lamu1./fu1 - lamu2./fu2;
  sig2 = lamu1./fu1 - lamu2./fu2;
  sigx = sig1 - sig2.^2./sig1;
  
  if (largescale)
    w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));
    h11pfun = @(z) -A(1./sigx.*At(z));
    [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;
    Adx = A(dx);
    Atdv = At(dv);
  else
    w1p = -(w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1)));
    H11p = A*(sparse(diag(1./sigx))*A');
    opts.POSDEF = true; opts.SYM = true;
    [dv,hcond] = linsolve(H11p, w1p, opts);
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;
    Adx = A*dx;
    Atdv = A'*dv;
  end
  
  du = (w2 - sig2.*dx)./sig1;
  
  dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
  dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;
  

  indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);
  s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
  indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);
  s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);

  suffdec = 0;
  backiter = 0;
  while (~suffdec)
    xp = x + s*dx;  up = u + s*du; 
    vp = v + s*dv;  Atvp = Atv + s*Atdv; 
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up;  fu2p = -xp - up;  
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    suffdec = (norm([rdp; rcp; rpp]) <= (1-alpha*s)*resnorm);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
      disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')
      xp = x;
      return
    end
  end

  x = xp;  u = up;
  v = vp;  Atv = Atvp; 
  lamu1 = lamu1p;  lamu2 = lamu2p;
  fu1 = fu1p;  fu2 = fu2p;

  sdg = -(fu1'*lamu1 + fu2'*lamu2);
  tau = mu*2*N/sdg;
  rpri = rpp;
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
  rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
  resnorm = norm([rdual; rcent; rpri]);
  done = (sdg < pdtol) | (pditer >= pdmaxiter);
end
