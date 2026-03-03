%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%????????ROI????????????voxel????????????????.
%??????????????????????????path??????????????????seed????mask
%????????????1*1*1????????????subdri??????????4D????.????????
%??path??????????????????????????????????monkey_prob_mask??.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

path = 'E:\cerebellar_MDTB';
path_data = 'E:\cerebellar_MDTB\FunImgARWglobalCFS';
path_result = 'E:\cerebellar_MDTB\Results\FC_FunImgARWglobalCFS_LN2Whole\';
path_seed ='E:\cerebellar_MDTB\LNperserved_lesion\';
mkdir(path_result);

subject = textread(strcat(path,filesep,'listApha.txt'),'%s');
%????????(????voxel)??????.
brain_mask = strcat(path,filesep,'Reslice_GreyMask.nii');
info = load_untouch_nii(brain_mask);
img = info.img;
[m n p] = size(img);
coordinates = zeros(1,1);
z = 1;
for i = 1:m
    for j = 1:n
        for k = 1:p
            if img(i,j,k) == 1
               coordinates(z,1) = i;
               coordinates(z,2) = j;
               coordinates(z,3) = k;
               z = z + 1;
            end
        end
    end
end
%save coordinate.mat coordinates
for sh = 1:length(subject);
    
seed_mask = strcat(path,filesep,'LNperserved_lesion',filesep,'LNPerserved_',subject{sh},'_LR.nii');
info = load_untouch_nii(seed_mask);
    img = info.img;
    [m n p] = size(img);
    coordinates_seed = zeros(1,1);
    z = 1;
    for ii = 1:m
        for j = 1:n
            for k = 1:p
                if img(ii,j,k) > 0
                   coordinates_seed(z,1) = ii;
                   coordinates_seed(z,2) = j;
                   coordinates_seed(z,3) = k;
                   z = z + 1;
                end
            end
        end
    end
    %save coordinates_WA2.mat coordinates_WA2

                                                               %%%%%%%%%%%%% revise the nki or blind name
        dirname = subject{sh};                                                              %%%%%%%%%%%%% revise the nki or blind name
        %subpath1 = strcat(path_nki,'\');
        subpath2= strcat(path_data,filesep,dirname);                                                  %%%%%%%%%%%%% revise the format
        subpath2dir1 = strcat(subpath2,filesep,'sFiltered_4DVolume.nii');
        info = load_untouch_nii(subpath2dir1);
        img_fMRI = info.img;
        [m n p q] = size(img_fMRI);
        %??????????????????voxel????????????.
        Time_course_target = zeros(length(coordinates),q);
        for wjj = 1:length(coordinates(:,1))
            Time_course_target(wjj,:) = img_fMRI(coordinates(wjj,1),coordinates(wjj,2),coordinates(wjj,3),:);
        end
    
        %??????????(????voxel)??????????.
        Time_courses = zeros(length(coordinates_seed(:,1)),q);
        for jjw = 1:length(coordinates_seed(:,1))
            Time_courses(jjw,:) = img_fMRI(coordinates_seed(jjw,1),coordinates_seed(jjw,2),coordinates_seed(jjw,3),:);
        end
        Time_courses_seed = mean(Time_courses);
        ws = length(Time_course_target);
        covariance_matrix = zeros(1,ws);
        for wh = 1:ws
            %covariance_seed_target = corrcoef(Time_courses_target{i},Time_course_seed(j,:));
            [r p] = corr(Time_courses_seed',(Time_course_target(wh,:))');
            covariance_matrix(1,wh) = r;
        end
        lovesh = find(isnan(covariance_matrix));
        covariance_matrix(lovesh) = 0;
        
        %save corrcoef.mat covariance_matrix
        %??????????Fisher????
        [m n] = size(covariance_matrix);
        covariance_matrix_Fisher = zeros(m,n);
        for jjsh = 1:n
            covariance_matrix_Fisher(1,jjsh) = 0.5 * log10((1 + covariance_matrix(1,jjsh))./(1 - covariance_matrix(1,jjsh)));
        end
        %????????????????????.
        %coor_mat = load('C:\Program Files\MATLAB\R2008a\work\functional_network\TimeCourses\coordinate.mat');
        %coordinates = coor_mat.coordinates;
        %template = strcat(path,'EPI_brain_mask_3mm.nii');
        info = load_untouch_nii(brain_mask);
        img = info.img;
        img(:,:,:) = 0;
        [m n] = size(covariance_matrix_Fisher);
        img = double(img);
        vector = covariance_matrix_Fisher;
        for wjjsh = 1:length(coordinates(:,1))
            img(coordinates(wjjsh,1),coordinates(wjjsh,2),coordinates(wjjsh,3)) = vector(wjjsh);
        end
        info.hdr.dime.datatype = 16;
        info.hdr.dime.bitpix = 32;
        info.img = img
        filename = strcat(path_result,subject{sh},'_LN_to_whole_brain');            %%%%%%%%%%%%% revise the nki or blind name                           
        save_untouch_nii(info,filename)
        disp(sh);

end


%%%% extrcat ROI mean value 
clear
clc

path = 'E:\cerebellar_MDTB';
pathgray = 'E:\cerebellar_MDTB\Results\FC_FunImgARWglobalCFS_LN2Whole';
path_seed ='E:\cerebellar_MDTB\';
list = textread(strcat('listApha.txt'),'%s');

for sh = 1:length(list)
    subseed = strcat(path,filesep,'rMDTB9mask.nii');

    info = load_untouch_nii(subseed);
    imgs = info.img;
    imgs = double(imgs);
    imgs(isnan(imgs)) = 0;
    findid = find(imgs>0);   
    file_1st = strcat(pathgray,filesep,list{sh},'_LN_to_whole_brain.nii');
    info = load_nii(file_1st);
    img1st = info.img;
    maxv1 = max(img1st(:));
    img1st = img1st./maxv1;
    tem_vec1 = img1st(findid);
    data2 = mean(tem_vec1);
    data3 = data2';
fid = fopen('meanvalue_rMDTB9mask.txt','a+');
fprintf(fid,'%15.7f\n',data3);
fclose(fid);
end
