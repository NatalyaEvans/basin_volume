%%

% This code is designed to load in bathymetry of Cariaco basin from GEBCO
% 2022 and estimate the volume between two depths. It uses a downloaded
% product from GEBCO_2022 Grid with a filename of
% 'gebco_2022_n45.1033_s40.8_w28.2_e41.9194.nc' as well as the accessory
% 'haversine' and its dependents. It first calculates the volume using the
% grid spacing on the GEBCO database, then interpolates at 1-2 orders of
% magnitude higher and re-integrates. Answers are printed at the end.

% To adjust the calculation, go down to the second section.

% Written 2023 04 11 by Talia Evans
% Adapted 2023 04 15 to work with Cariaco basin

%% load data

clear all
% ncdisp('gebco_2022_n11.2_s10.1_w-66.25_e-64.0.nc')
lat=ncread('gebco_2022_n11.2_s10.1_w-66.25_e-64.0.nc','lat'); % lat in degree N
long=ncread('gebco_2022_n11.2_s10.1_w-66.25_e-64.0.nc','lon'); % lat in degree E
z=ncread('gebco_2022_n11.2_s10.1_w-66.25_e-64.0.nc','elevation'); % elevation in m. Note that depth is negative

% plot data real quick
figure(1)
surf(long,lat,z','LineStyle','none')
colorbar
xlabel(['Longitude/' char(176) 'E'])
ylabel(['Latitude/' char(176) 'N'])
zlabel('Elevation/m')
title('Downloaded Cariaco basin bathymetry')

% contourf(long,lat,z',[-700:10:0],'LineStyle','none')
% contourf(long,lat,z','LineStyle','none')

% colorbar

%% select bounds 

% DEPTHS ARE NEGATIVE

% zbot=min(z,[],'all'); % deepest point
% ztop=double(-900); 

% zbot=double(-900); 
% ztop=double(-550);

zbot=double(-550); 
ztop=double(-350);

% long_cut_check=1; % if we have to cut the longitude into a section, set this to one
long_cut_check=0; % if we have to cut the longitude into a section, set this to one
long_cut=-65.1;
long_side=1; % 0 for W, 1 for E basin

% shal_clean=0; % set this to 1 if integrating to 350 m to clean the NE corner
shal_clean=1; % set this to 1 if integrating to 350 m to clean the NE corner

high_res=1; % if high res=1, then a high res integration is calculated. If not one, this section is skipped
% high_res=0; % if high res=1, then a high res integration is calculated. If not one, this section is skipped
res=200000-1; % resolution (number of bins) for high res integration
interptype='spline'; % interpolation type for high res integration. Options are linear, pchip, spline, and nearest

%% resize data, filter, integrate

ind= z >= zbot & z <= ztop; % identify relevant elevation values
ind_lat=[sum(ind,1)~=0]'; % index for y range
ind_long=[sum(ind,2)~=0]; % index for x range

z2=z(ind_long,ind_lat); % subset elevation
lat2=lat(ind_lat); % subset lat
long2=long(ind_long); % subset long

if long_cut_check==1
    ind= long2<long_cut;
    if long_side==0
        long2(ind==1)=[];
        z2(ind==1,:)=[];
    end
    if long_side==1
        long2(ind==0)=[];
        z2(ind==0,:)=[];
    end
end      

m=length(long2);
n=length(lat2);

if shal_clean==1
    longz=repmat(long2,[1,n]);
    latz=repmat(lat2,[1,m])';
    ind= [latz > [-0.2.*longz-2.1]] & z2<ztop;
    z2(ind)=-299;
    
    ind= [latz > 11] & z2<ztop;
    z2(ind)=-299;
end



lat_m=ones(size(n)).*NaN; % create grid to fill in
long_m=ones(size(m)).*NaN; % create grid to fill in
long_ref=min(long2); % define reference for distance
lat_ref=min(lat2); % define reference for distance

for i=1:n
    lat_m(i)=1000*haversine([long_ref,lat2(i)],[long_ref,lat_ref]); % distance for the lat in m
end

for j=1:m
    long_m(j)=1000*haversine([long2(j),lat_ref],[long_ref,lat_ref]); % distance for the longitude in m
end

% mask data integral

z3=double(z2); % duplicate data
% z3(z3>ztop)=NaN; % mask data that is too shallow but in the range
z3(z3>ztop)=ztop; % mask data that is too shallow but in the range
z3(z3<zbot)=zbot; % mask data that is too shallow but in the range


figure(2)
surf(long2,lat2,z3','LineStyle','none');
% contourf(long2,lat2,z3','LineStyle','none');
colorbar
xlabel(['Longitude/' char(176) 'E'])
ylabel(['Latitude/' char(176) 'N'])
zlabel('Elevation/m')
title('Selected Cariaco basin bathymetry')

% perform integral across lat
int1=ones(1,m).*NaN;
for j=1:m
    z4=z3(j,isnan(z3(j,:))==0); % make sure data is clean
    if length(z4)>=2
        lat_m2=lat_m(isnan(z3(j,:))==0);
        int1(j)=trapz(lat_m2,z4-ztop); % integrate across latitude, m^2, zeroing the reference surface
    end
end

long_m2=long_m(isnan(int1)==0); % just a check for NaN
int2=int1(isnan(int1)==0); 
int_m=-trapz(long_m2,int2); % integrate across longitude, m^3

display('Upper and lower bounds for integration')
[ztop zbot]
display('Integrated volume/km^3 at grid size resolution')
int_km=int_m/(1000^3) % int in km^3


%% high res
% Try a high resolution to see if estimates change

if high_res==1
    
    lat_mg=[0:max(lat_m)/res:max(lat_m)];
    for j=1:m
        z_interp=interp1(lat_m,double(z2(j,:)),lat_mg,interptype);
        z3=z_interp;
        %     z3(z_interp>ztop)=NaN; % mask value
        %     z4=z3(isnan(z3)==0);
        %     lat_mg2=lat_mg(isnan(z3)==0);
        lat_mg2=lat_mg;
        z3(z_interp>ztop)=ztop; % zero value, which keeps horizontal resolution
        z3(z_interp<zbot)=zbot; % zero value, which keeps horizontal resolution
        
        z4=z3;
        if length(z4)>=2
            inth1(j)=trapz(lat_mg2,z4-ztop); % integrate across latitude, m^2
        end
    end
    
    long_m2=long_m(isnan(inth1)==0);
    inth2=inth1(isnan(inth1)==0);
    
    long_mg2=[0:max(long_m)/res:max(long_m)];
    inth2_interp=interp1(long_m,inth2,long_mg2,interptype);
    
    inth_m=-trapz(long_mg2,inth2_interp); % integrate across longitude, m^3
    
    display('Integrated volume/km^3 at higher resolution')
    inth_km=inth_m/(1000^3) % int in km^3
    
    display('Resolution increase in directions x and y')
    [res./length(long_m) res./length(lat_m)]
    
    display('Interpolation technique')
    interptype
end
