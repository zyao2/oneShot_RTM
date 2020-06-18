function seis=forward(v,nbc,dx,nt,dt,s,sx,sz,gx,gz)
seis=zeros(nt,numel(gx));[nz,nx]=size(v);
ng=numel(gx);
c1 = -205.0/72.0; c2 = 8.0/5.0; c3 = -1.0/5.0; c4 = 8.0/315.0; c5 = -1.0/560.0;

% setup ABC and temperary variables
v=padvel(v,nbc);
abc=Coef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*c1*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);
[isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc);
p1=zeros(size(v)); p0=zeros(size(v));

% Time Looping
for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (c2*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +c3*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +c4*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +c5*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
    % dipole source
    %p(isz-2,isx) = p(isz-2,isx) - beta_dt(isz-2,isx) *wavelet(it);
    if isFS
        p(nbc+1,:)=0.0;
        p(nbc:-1:nbc-3,:) = - p(nbc+2:nbc+5,:);
    end
    for ig=1:ng
        seis(it,ig)=p(igz(ig),igx(ig));
    end
    if nargout>1 % save BC
        [bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it)]=save_boundary(p,nz,nx,nbc);
    end
    p0=p1;
    p1=p;
end
