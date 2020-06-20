%Single shot RTM;
%velocity model
clear all
nz=81;nx=201;
vel=zeros(nz,nx);
vel(1:30,:)=1000;
vel(31:60,:)=1200;
vel(61:end,:)=1500;

%source and receiver location;
dx=5;dz=5;
sx=100*dx;sz=0;
x = (0:nx-1)*dx; z = (0:nz-1)*dx;
recx=(0:2:(nx-1))*dx; recz=zeros(size(recx));
nrec=length(recx);
%FD parameters;
nbc=20; nt=2001; dt=0.0005; t=(0:nt-1)*dt;

%source wavelet;
freq=25; s=ricker(freq,dt);
% seismic data by  modeling with true velocity model;
%seis=forward(vel,nbc,dx,nt,dt,s,sx,sz,recx,recz);

load seis

%Smooth the true veolicyt to get the migration velocity model;

[vel_ss,refl_ss]=vel_smooth(vel,3,3,1);

%Plot the velocity and seismic data;

figure(1);set(gcf,'position',[0 0 1000 400]);subplot(231);imagesc(x,z,vel);colorbar;
xlabel('X (m)'); ylabel('Z (m)'); title('velocity');
figure(1); subplot(232);imagesc(1:nrec,t,seis);colormap(gray);
title('Seismic Profile');ylabel('Time (s)');xlabel('rec #');caxis([-0.25 0.25]);

%Run the rtm code, watch the movie of forward propagation field, reconstruction forward field, backward propagation field and accumulated rtm image;

tic;
[img,illum]=rtm(seis,vel_ss,nbc,dx,nt,dt,s,sx,sz,recx,recz);
toc;

%Plot the illumination compensation, and see the impact.

figure(2); set(gcf,'position',[0 0 400 600]);colormap(gray);
subplot(311);imagesc(x,z,img);caxis([-10 10]);
xlabel('X (m)'); ylabel('Z (m)'); title('rtm image');
figure(2);subplot(312);imagesc(x,z,illum);caxis([-100 100]);
xlabel('X (m)'); ylabel('Z (m)'); title('illumination compensation');
figure(2);subplot(313);imagesc(x,z,img./illum);caxis([-1 1]);
xlabel('X (m)'); ylabel('Z (m)'); title('rtm image after compensation');


