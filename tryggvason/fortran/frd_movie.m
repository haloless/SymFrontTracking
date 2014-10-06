%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A Simple program to plot the output from FTC2D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fileName='ttt'
% fileName='ex1'
fileName='ex2'
numberOfFiles=50
pl_frd=1            % plot front 1=yes
pl_vel=1           % plot velocity 1=yes
pl_rho=0            % plot density 1=yes (only if density fields are printed)

xl=2.0; yl=2.0; %Size of domain

c1=0.25*xl; %arbitrary constant

for fileNumber=1:numberOfFiles

  if(pl_frd==1)
  	inf1=sprintf('/frd%03d',fileNumber);
	infile=sprintf([fileName,inf1])
	[fid, message]=fopen(infile,'r')	
	nn=fscanf(fid,'%d',[4]);
	np=nn(2)
	ne=nn(3);
	n1=fscanf(fid,'%d',[2]);
	pp=fliplr(rot90(fscanf(fid,'%f',[2,np]),3));
	ll=fliplr(rot90(fscanf(fid,'%d',[2,ne]),3));	
	fclose(fid);

% put bubbles back into periodic box 
   for i=1:np,xx=pp(i,1);pp(i,1)=xx-xl*floor(xx/xl);end
   for i=1:np,yy=pp(i,2);pp(i,2)=yy-yl*floor(yy/yl);end
	
	plot(0,0,'wo');hold on
	for i=1:ne	
	   if( (pp(ll(i,1),1)-pp(-ll(i,2),1))^2+(pp(ll(i,1),2)-pp(-ll(i,2),2))^2 < c1)
	       plot([pp(ll(i,1),1),pp(-ll(i,2),1)],...
	             [pp(ll(i,1),2),pp(-ll(i,2),2)],'y','LineWidth',2);
	    end
     end
  end

  if(pl_vel==1)  
	inf1=sprintf('/uvl%03d',fileNumber);
	infile=sprintf([fileName,inf1])
	[fid, message]=fopen(infile,'r');
	line1=fgetl(fid);	
	nn=fscanf(fid,'%d',[2]);
	nx=nn(1);
	ny=nn(2);
	xc=fscanf(fid,'%f',[nn(1)]);
	yc=fscanf(fid,'%f',[nn(2)]);
	u=fliplr(rot90(fscanf(fid,'%f',[nx,ny]),3));
	fclose(fid);

	inf1=sprintf('/vvl%03d',fileNumber);
	infile=sprintf([fileName,inf1])
	[fid, message]=fopen(infile,'r');
	line1=fgetl(fid);	
	nn2=fscanf(fid,'%d',[2]);
	nx2=nn2(1);
	ny2=nn2(2);
	xc2=fscanf(fid,'%f',[nn2(1)]);
	yc2=fscanf(fid,'%f',[nn2(2)]);
	v=fliplr(rot90(fscanf(fid,'%f',[nx2,ny2]),3));
	fclose(fid);	

    uav(1:ny-1,1:nx-1)=0.5*(u(1:ny-1,1:nx-1)+u(2:ny,1:nx-1));
    vav(1:ny-1,1:nx-1)=0.5*(v(1:ny-1,1:nx-1)+v(1:ny-1,2:nx));

	plot([0 xc(nx) xc(nx) 0 0],[0 0 yc(ny2) yc(ny2) 0]);
	axis 'equal', axis([0 xc(nx) 0 yc(ny2)])

    xt(1:nx-1)=xc(1:nx-1);yt(1:ny-1)=yc(1:ny-1);	
	quiver(xt,yt,uav,vav);
  end

  if(pl_rho==1)  
	inf1=sprintf('/rho%03d',fileNumber);
	infile=sprintf([fileName,inf1])
	[fid, message]=fopen(infile,'r');
	line1=fgetl(fid);	
	nn3=fscanf(fid,'%d',[2]);
	nx3=nn3(1);
	ny3=nn3(2);
	xc3=fscanf(fid,'%f',[nn3(1)]);
	yc3=fscanf(fid,'%f',[nn3(2)]);
	r=fliplr(rot90(fscanf(fid,'%f',[nx3,ny3]),3));
	fclose(fid);	
%    contour(xc3,yt,r)
	mesh(r)
  end

%  eval(['print -dpict ','xxx',num2str(fileNumber)]);    %to make a movie
  hold off
  pause;
end