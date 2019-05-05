% Test demo
    % A test demo of detecting clustered microcalcifications on DBT volumes.
    % Author: Zhang Fan, E-mail:zf2016@mail.ustc.edu.cn
    % https://github.com/zhangfan2018
    % Fan Zhang,et al. "Multi-domain features for reducing false positives 
    % in automated detection of clustered microcalcifications in digital 
    % breast tomosynthesis." Medical Physics, 2019.
% Last modified 2019.02.27

close all; clear; clc;
warning off;
functiondir=which('TestDemo.m'); functiondir=functiondir(1:end-length('TestDemo.m'));
addpath([functiondir '/Library'],[functiondir '/breastImages'],[functiondir '/annoationFiles']);

% Reading orignal DBT volume and the corresponding annotation file.
[imageData,~]=nrrdread('1605120839_L_CC.nrrd');
labelInf=labelInfRead('1605120839_L_CC.xml');
lesion_center=str2num(labelInf.(['zCenter',num2str(1)]));

% DBT volume is cropped based on the middle slice.
[croped_image,breast_region,cropped_row,breast_boundary,dis_nipple]=CropImage(imageData);
% Plotting the breast boundary.
[~,N]=size(breast_region);
for i=1:1
    figure;imshow(croped_image(:,:,lesion_center),[]);hold on;
    for j=1:N-1
        plot([breast_region(lesion_center,j),breast_region(lesion_center,j+1)],...
            [j,j+1],'b','LineWidth',2);
    end
end

% Enhancementing the contrast-to-noise ratio of microcalcification.
cnr_image=EnhancementCnr(croped_image,breast_region);
figure;imshow(cnr_image(:,:,lesion_center),[]);

% Detecting the individual microcalcification, including its location and
% boundary.
[object_seeds,objects_number]=DetectObject(cnr_image,breast_region);

% Plotting the location of detected objects.
object_center_tmp=regionprops3(object_seeds,'Centroid');
object_center=round(cat(1, object_center_tmp.Centroid));
object_center(:,[2,1])=object_center(:,[1,2]);
PlotLocMap(croped_image,object_center,labelInf,cropped_row);
assert(objects_number~=0, 'There is not the detected microcalcification object');

% Segmenting the individual microcalcification is to obtain the accurate boundary.
% [segmented_mask,object_center,ind_3D_features]=SegmentObject(cnr_image,object_seeds,breast_region,breast_boundary);
[segmented_mask,object_center,ind_features]=activeSegment(cnr_image,object_seeds,breast_region,breast_boundary);

% Displaying the boundary of detected objects.
figure;imshow(croped_image(:,:,lesion_center),[200 750]);hold on;
visboundaries(object_seeds(:,:,lesion_center),'Color','b','LineWidth',1); 
% Displaying the boundary of segmented objects.
figure;imshow(croped_image(:,:,lesion_center),[200 750]);hold on;
visboundaries(segmented_mask(:,:,lesion_center),'Color','b','LineWidth',1); 

% Plotting the location of segmented objects.
PlotLocMap(croped_image,object_center,labelInf,cropped_row);

% Display the 3D model for the segmented objects.
figure;
p = patch(isosurface(double(segmented_mask)));
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 27/128]);
camlight; 
lighting phong

% Danamic clustering
% [cluster_result_center,cluster_result_scope,cluster_result_temp,cluster_features,focus_slice_z] = ...
%     RegionCluster(croped_image,segmented_mask,object_center,ind_3D_features,breast_region,breast_boundary,dis_nipple);
[cluster_result_center,cluster_result_scope,cluster_result_temp,cluster_features] = ...
    dynamicCluster(croped_image,segmented_mask,object_center,ind_features,breast_region,breast_boundary,dis_nipple);

% % Feature nomalization.
% [m_temp,n_temp]=size(cluster_features);
% features=cluster_features(:,10:331);
% load ps;
% train_scale=mapminmax('apply',features',ps);
% train_scale=train_scale';
% 
% % Data transformation.
% m_flag=0;
% load featureslda;
% feature_vector=features;
% clear features m_data_temp;
% for i=1:332
%     for j=1:26
%     if i==feature_vector(j)
%         if m_flag==0
%             m_data_temp=train_scale(:,i);
%             m_flag=1;
%         else
%             m_data_temp=[m_data_temp train_scale(:,i)];
%         end
%     end
%     end
% end
% 
% % FPs reduction by LDC classifier.
% load modelLda;
% [~,dec_values_temp] = predict(MdlLinear,m_data_temp);
% m_cnt=0;
% clear cluster_result_center_temp cluster_result_scope_temp index_temp;
% for k=1:m_temp
%     if dec_values_temp(k,2)>0
%         m_cnt=m_cnt+1;
%         cluster_result_center_temp(m_cnt,:)=cluster_result_center(k,:);
%         cluster_result_scope_temp(m_cnt,:)=cluster_result_scope(k,:);
%         index_temp(m_cnt,1)=k;
%     end
% end

% Display the detection results.
if cluster_result_center~=0
    [M,~]=size(cluster_result_center);
    [N,~]=size(cluster_result_temp);
    figure;
    imshow(croped_image(:,:,lesion_center),[200 750]);
    hold on;
    for i=1:M
        for j=1:N
            if cluster_result_temp(j,4)==i
                plot(cluster_result_temp(j,2),cluster_result_temp(j,1),'b^','LineWidth',1,'markersize',3);
                hold on;
            end
        end
        x=cluster_result_scope(i,3)-5;
        y=cluster_result_scope(i,1)-5;
        width=cluster_result_scope(i,4)-cluster_result_scope(i,3)+14;
        height=cluster_result_scope(i,2)-cluster_result_scope(i,1)+14; 
        rectangle('Position',[x,y,width,height],'LineStyle','--','Curvature',[1 1],'LineWidth',1,'EdgeColor','b')
        a=double(x+width/2-20);
        b=double(y+height+10);
        text('Position',[a,b],'string','score=0.93','color','red','FontSize',5);
    end
end

