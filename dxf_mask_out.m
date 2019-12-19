function [  ] = dxf_mask_out( crystal_mask ,file_name,L_HW,L_y,L_x,dy)
%% Output DXF file
% 1 is poled region
PP1=zeros([2,1].*size(crystal_mask));
PP2=PP1;
% PP1(1:2:end)=~crystal_mask;
% PP2(2:2:end)=crystal_mask;
PP1(1:2:end)=(crystal_mask==0);
PP2(2:2:end)=(crystal_mask==1);
Pixels1 = bwboundaries(PP1);
Pixels2 = bwboundaries(PP2);

m_max=cell2mat(cellfun(@max,Pixels1,'UniformOutput', false));
m_min=cell2mat(cellfun(@min,Pixels1,'UniformOutput', false)); 
m_mix_1=[m_min(:,1),m_max(:,2)]; % x-min,y-max
m_mix_2=[m_max(:,1),m_min(:,2)];% x-max,y-min
m1=[m_min,m_mix_1,m_max,m_mix_2];
m1=reshape(m1,[size(m1,1),2,4]);
m_max=cell2mat(cellfun(@max,Pixels2,'UniformOutput', false));
m_min=cell2mat(cellfun(@min,Pixels2,'UniformOutput', false)); 
m_mix_1=[m_min(:,1),m_max(:,2)]; % x-min,y-max
m_mix_2=[m_max(:,1),m_min(:,2)];% x-max,y-min
m2=[m_min,m_mix_1,m_max,m_mix_2];
m2=reshape(m2,[size(m2,1),2,4]);
m=cat(1,m1,m2);
% figure;scatter(reshape(m(:,1,:),[],1),reshape(m(:,2,:),[],1))
% shifted center location
shift_center=round(size(PP1)/2);

%opening the final DXF file
for mult_factor=[1,0.5,0.25]
    filename=[file_name,num2str(100*mult_factor),'.dxf'];
    [fid,~]=DXF_start(filename,1);
    number_of_polygons=size(m,1);
    for j=1:number_of_polygons
        y = squeeze((m(j,1,:)-shift_center(1))*L_HW)+[-1;-1;1;1]*L_HW/2;
        y(3:4) = (y(3:4)-y(1:2))*mult_factor+y(1:2); % Scale duty cycle
        y = [y;y(1)];
        x = squeeze((m(j,2,:)-shift_center(2))*dy)+[-1;1;1;-1]*dy/2;
        x= [x;x(1)];
        DXF_poly(fid,x,y,length(x),7,1);
        disp(['SG ' num2str(j/number_of_polygons*100),'%']);
    end
    % Create the frame
%     frame_thickness=0.1*1e3;
%     x0=0; y0=0;
%     mask_gap=0*0.1e3;%0.2e3;
%     Rinx=L_y/2+mask_gap; Routx=Rinx+frame_thickness;
%     Riny=L_x/2+mask_gap; Routy=Riny+frame_thickness;
%     points=poly_donut_closed_advanced(x0,y0,Rinx,Routx,Riny,Routy);
%     for ii=1:4
%         y=points(ii).y;
%         x=points(ii).x;
%         S=length(x);
%         DXF_poly(fid,x,y,S,7,1);
%     end
    DXF_end(fid);
    fclose('all');
end
end

