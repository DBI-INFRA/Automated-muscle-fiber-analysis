
function muscleFibreAnalysis_for_review(aImarisApplicationID, demo_data_path)
%
% function muscleFibreAnalysis(aImarisApplicationID)
%
%   Input:
%   aImarisApplicationID        Imaris server object, representing the
%                               currently active scene in Imaris. 
%                               The ID is provided by Imaris when calling 
%                               the plugin.
%
% muscleFibreAnalysis is a Matlab plugin for the Imaris software package 
%   (https://imaris.oxinst.com/open/). It queries the current scene for an
%   Imaris surface object (muscle fiber) and a set of Imaris spots objects
%   (myonuclei), imports the data to Matlab and performs a series of 
%   measurements:
%
%   1.) Inside nuclei NN [um]
%       Nearest neighbor (NN) distance between all myonuclei labeled as 
%       being located inside the fiber
%   2.) Outside nuclei NN [um] 
%       NN distance between all nuclei labeled as being located outside 
%       the fiber
%   3.) Outside to Inside nuclei NN [um]
%       Inside nuclei NN distance for all outside nuclei
%   4.) Outside to fiber surface distance [um]
%       Shortest distance to the fiber surface for all outside nuclei
%   5.) Myonuclei domain surface area (MNDS) [um^2]
%       The area of the fiber surface, that is part of a nucleus' domain.
%       Areas touching the data domain boundary are discarded.
%   6.) Myonuclei domain volume (MNDV) [um^3]
%       The volume of the fiber, that is part of a nucleus' domain. Volumes
%       touching the data domain boundary are discarded.
%   7.) The number of outside and inside nuclei per unit slice [count],
%       measured with respect to the longitudinal axis of the fiber, 
%       starting from the tip of the fiber. The slice thickness can be set 
%       by the user. Data slices touching the data domain boundary are 
%       discarded. 
%   8.) Inside nuclei density [count/um^3] 
%       Number of inside nuclei per slice volume.
%
%   All analysis results are provided as Matlab figures and exported as PNG 
%   image files and Excel sheets. In addition, the following figures are 
%   generated and exported:
%   a.) The fiber surface and the inside and outside nuclei
%   b.) The MNDS domains in different colors
%   c.) The MNDS domains, semi-transparent, and the MNDV sizes represented
%       as spheres around the nuclei
%   d.) A 3D image stack with the MNDVs in different greyvalues
%
%   IMPORTANT: Nuclei inside the fiber should be labeled 'Class A' and 
%   nuclei outside the fiber should be labeled 'Class B' in Imaris, as the
%   plugin cannot distinguish between the two otherwise.
%
%
% Version: 1.0
% Author: Martin Baiker-Soerensen; DBI-IACF; Copenhagen University; 03/23
% Update: Martin Baiker-Soerensen; DBI-IACF; Copenhagen University; 04/23
% Status: TESTING (Matlab version: 2022b, Imaris version: 10.0.0)
%
% Note that the code related to the Imaris bridge communication was in part
% developed by Aaron Ponti in 2012, as part of an Imaris-Matlab bridge
% workshop: http://www.scs2.net/next/files/courses/iic/ImarisXTCourse.pdf
%
% Note that the code for aligning the sample to the longitudinal axis is 
% based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z 
% direction based on PCA, MATLAB Central File Exchange. April 21, 2023.'
%
% Online documentation for the project:
% https://alumni-my.sharepoint.com/:w:/g/personal/tmj555_ku_dk/EdZmaiEQ_mhMjfABJxJKGfABR6Zz-X5M_3gOBgFgEvpOiw?e=nbzaXn
% ---------------------------------------------------------
%
%   Imaris header (compulsory for the plugin to be added to the menu):
%
%    <CustomTools> 
%      <Menu> 
%        <Submenu name="Matlab-Imaris Bridge"> 
%          <Item name="Muscle Fiber Analysis" icon="Matlab"> 
%            <Command>MatlabXT::muscleFibreAnalysis(%i)</Command> 
%          </Item> 
%        </Submenu>
%      </Menu> 
%      <SurpassTab> 
%        <SurpassComponent name="bpMatlabImaris"> 
%          <Item name="Muscle Fiber Analysis" icon="Matlab"> 
%            <Command>MatlabXT::muscleFibreAnalysis(%i)</Command> 
%          </Item> 
%        </SurpassComponent> 
%      </SurpassTab> 
%    </CustomTools> 
%



%% For debugging only
show_interm_figures = 0;
debugging           = 0;



%% If used for code review, the input argument is the path to the folder with
%  the mat file including test data
if isempty(aImarisApplicationID)
    load(fullfile(demo_data_path, 'demo_data.mat'));
else
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up connection between Imaris and MATLAB and get data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % The Imaris server object should be passed to the plugin by Imaris. If
    % that didn't work somehow, look for the object
    if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
        vImarisApplication = aImarisApplicationID;
    else
        % Connect to Imaris interface
        javaaddpath ImarisLib.jar
        vImarisLib = ImarisLib;
        if ischar(aImarisApplicationID)
            aImarisApplicationID = round(str2double(aImarisApplicationID));
        end
        vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
    end
    
    
    % Get the Surfaces (ISurface objects in Imaris), volumes (IVolume objects 
    % in Imaris) and spots (ISpots objects in Imaris)
    surfaces    = {};   spots       = {};   volumes = {}; 
    nSurfaces   = 0;    nVolumes    = 0;    nSpots  = 0;
    
    % Get current scene objects from Imaris and find out what we got
    nChildren = vImarisApplication.GetSurpassScene().GetNumberOfChildren();
    for i = 0 : (nChildren - 1)
        child = vImarisApplication.GetSurpassScene.GetChild( i );
        if vImarisApplication.GetFactory().IsSurfaces(child)
            nSurfaces = nSurfaces + 1;
            % We must cast the child to a Surfaces object!
            surfaces{nSurfaces} = ...
                vImarisApplication.GetFactory().ToSurfaces(child);
        end
        if vImarisApplication.GetFactory().IsVolume(child)
            nVolumes = nVolumes + 1;
            % We must cast the child to a Volume object!
            volumes{nVolumes} = ...
                vImarisApplication.GetFactory().ToVolume(child);
        end 
        if vImarisApplication.GetFactory().IsSpots(child)
            nSpots = nSpots + 1;
            % We must cast the child to a Spots object!        
            spots{nSpots} = ...
                vImarisApplication.GetFactory().ToSpots(child);
        end
    end
    
    
    
    %% The spots represent the nuclei. Get their Centers of Gravity (CoGs)
    CoGs = spots{1}.GetPositionsXYZ;
    
    % Get all spot class labels
    % IMPORTANT: Nuclei inside the fiber should be labeled 'Class A' and 
    %   nuclei outside the fiber should be labeled 'Class B' in Imaris, as the
    %   plugin cannot distinguish between the two otherwise. 
    tmp          = spots{1}.GetLabels;
    class_labels = cell(size(CoGs, 1), 1);
    CoGs_inside  = [];  % Class A
    CoGs_outside = [];  % Class B
    
    for CoG_cnt = 1:size(CoGs, 1)
        class_labels{CoG_cnt} = tmp(CoG_cnt).mLabelValue;
    
        if strcmp(tmp(CoG_cnt).mLabelValue, 'Class A')
            CoGs_inside  = [CoGs_inside; CoGs(CoG_cnt,:)];
        elseif strcmp(tmp(CoG_cnt).mLabelValue, 'Class B')
            CoGs_outside = [CoGs_outside; CoGs(CoG_cnt,:)];
        else
            error('Class label not know. Has to be ''Class A'' for spots inside and ''Class B'' for spots outside the fiber.' );
        end
    end
    
    
    
    %% Extract actual surface from Imaris ISurface object
    cur_surface     = surfaces{1};                          % ISurface
    % Get surface data (note that we need to use Java indexing now...)
    cur_surf_data   = cur_surface.GetSurfaceData(0);        % IDataSet
    % Get surface dimensions
    surface_layout  = cur_surface.GetSurfaceDataLayout(0);
    % Generate mask from data
    surf_mask       = cur_surface.GetMask(surface_layout.mExtendMinX, ...
                                          surface_layout.mExtendMinY, ...
                                          surface_layout.mExtendMinZ, ...
                                          surface_layout.mExtendMaxX, ...
                                          surface_layout.mExtendMaxY, ...
                                          surface_layout.mExtendMaxZ, ...
                                          surface_layout.mSizeX, ...
                                          surface_layout.mSizeY, ...
                                          surface_layout.mSizeZ, ...
                                          0);               % IDataSet
    
    % Get the actual data out of the IDataSet object... 
    % ... first get the dataset class ...
    switch char(cur_surf_data.GetType())
        case 'eTypeUInt8', datatype = 'uint8';
        case 'eTypeUInt16', datatype = 'uint16';
        case 'eTypeFloat', datatype = 'single';
        otherwise, error('Bad value for iDataSet.GetType()');
    end
    
    % ... allocate memory ...
    surf_mask_stack  = zeros([surf_mask.GetSizeX(), surf_mask.GetSizeY(), ...
                              surf_mask.GetSizeZ()], datatype);
    
    % ... get the data ...
    switch char(cur_surf_data.GetType())
        case 'eTypeUInt8'
            arr                 = surf_mask.GetDataVolumeAs1DArrayBytes(0, 0);
            surf_mask_stack(:)  = typecast(arr, 'uint8');
        case 'eTypeUInt16'
            arr = surf_mask.GetDataVolumeAs1DArrayShorts(0, 0);
            surf_mask_stack(:)  = typecast(arr, 'uint16');
            surf_mask_stack(:)  = uint8(double(surf_mask_stack(:))./double(max(surf_mask_stack(:))).*255);
            surf_mask_stack     = uint8(surf_mask_stack);            
        case 'eTypeFloat'
            surf_mask_stack(:)  = surf_mask.GetDataVolumeAs1DArrayFloats(0, 0);
        otherwise
            error('Bad value for type');
    end
    
    % ... and flip x and y dimensions (Imaris (Java) vs. Matlab...)
    surf_mask_stack = permute(surf_mask_stack, [2 1 3]);
    
    clear arr;
    
    if debugging
        figure('Name', 'Data loaded');
        pause(.1);
    else
        disp('Data successfully imported ...');
    end
    
    
    
    %% Determine the voxel size
    % ImarisXT does not offer to get the voxel size... so we have to calculate it
    % 'cur_surf_data' should be an iDataSet
    voxel_size_x = (cur_surf_data.GetExtendMaxY() - cur_surf_data.GetExtendMinY()) / ...
                   (cur_surf_data.GetSizeY() - 1);
    voxel_size_y = (cur_surf_data.GetExtendMaxX() - cur_surf_data.GetExtendMinX()) / ...
                   (cur_surf_data.GetSizeX() - 1);
    voxel_size_z = (cur_surf_data.GetExtendMaxZ() - cur_surf_data.GetExtendMinZ()) / ...
                   (cur_surf_data.GetSizeZ() -1);
    
    voxel_size   = [voxel_size_x voxel_size_y voxel_size_z];        % in um
    voxel_volume = voxel_size_x * voxel_size_y * voxel_size_z;      % in um^3

    % Adjust CoG values to current reference frame
    % A point at (0|0) corresponds after isosurface extraction to a point at (0.5|0.5)
    CoGs_inside = [ CoGs_inside(:,1)  - cur_surf_data.GetExtendMinX() + 0.5 ...
                    CoGs_inside(:,2)  - cur_surf_data.GetExtendMinY() + 0.5 ...
                    CoGs_inside(:,3)  - cur_surf_data.GetExtendMinZ() + 0.5];
    CoGs_outside = [CoGs_outside(:,1) - cur_surf_data.GetExtendMinX() + 0.5 ...
                    CoGs_outside(:,2) - cur_surf_data.GetExtendMinY() + 0.5 ...
                    CoGs_outside(:,3) - cur_surf_data.GetExtendMinZ() + 0.5];

    save(fullfile(demo_data_path, 'demo_data.mat'));
    return;
end



%% Ask user to provide a unique data ID, the slice thickness, a subsampling
%  factor, a path to a folder to store results and the name of the Excel
%  sheet where all data should be stored
slice_thickness_default         = '15';   % in um
subsampling_factor_default      = '1';
results_saving_folder_default   = 'F:\Christian_H\Fiber_analysis\Results';
results_file_default            = 'Christian_MTJ_MND_project_results.xlsx';

answer  = inputdlg({'Unique ID', 'Slice thickness (um)', 'Subsampling', 'Results folder', 'Results Excel filename'}, ...
                    'Specify parameters', [1 75], {'', slice_thickness_default, subsampling_factor_default, ...
                    results_saving_folder_default, results_file_default});

% If no ID was provided, an error occurs
if isempty(answer{1})
    f = errordlg('Please provide valide data ID!', ...
                 'Data ID missing...' ); 
    return;
else
    data_ID            = answer(1);
end

slice_thickness         = str2double(answer{2});
subsampling_factor      = str2double(answer{3});
results_saving_folder   = answer{4};
results_file            = answer{5};
excel_file_path         = fullfile(results_saving_folder, results_file);



%% Extract fiber surface
% If needed (typically not), subsample the volume, and then artificially 
% increase the volume to avoid holes in the areas touching the borders when
% extracting the surface using isosurface()
surf_mask_stack           = logical(surf_mask_stack);
surf_mask_stack           = surf_mask_stack(1:subsampling_factor:end, ...
                                            1:subsampling_factor:end, ...
                                            1:subsampling_factor:end); 

voxel_volume              = voxel_volume * subsampling_factor^3;

surf_mask_stack_enlarged = zeros(size(surf_mask_stack) + 2, 'logical');
surf_mask_stack_enlarged(2:end-1, 2:end-1, 2:end-1) = surf_mask_stack; 

% Extract the fiber surface
[x, y, z]   = meshgrid(1:size(surf_mask_stack_enlarged, 2), ...
                       1:size(surf_mask_stack_enlarged, 1), ...
                       1:size(surf_mask_stack_enlarged, 3)); 

fv = isosurface(x,y,z,surf_mask_stack_enlarged, 0.5);

% Correct for the enlarged volume size
fv.vertices = fv.vertices - 1.5;
fv_orig     = fv.vertices;

fv.vertices(:,1) = fv.vertices(:,1).*voxel_size(1).*subsampling_factor;
fv.vertices(:,2) = fv.vertices(:,2).*voxel_size(2).*subsampling_factor;
fv.vertices(:,3) = fv.vertices(:,3).*voxel_size(3).*subsampling_factor;

if debugging
    figure('Name', 'Surface extracted');
    pause(.1);
else
    disp('Surface extracted ...');
end

if show_interm_figures
    figure;
    colorfield = ones(length(fv.vertices(:,1)), 1);
    facecolor = [1 0.8 0.4];
    edgecolor = 'black';
    trimesh(fv.faces, fv.vertices(:,1), fv.vertices(:,2), ...
                      fv.vertices(:,3), ...
                      colorfield, 'FaceColor', facecolor, ...
                     'EdgeColor', edgecolor, 'FaceAlpha', 0.5);
    daspect([1 1 1])
    view(3); axis tight
    camlight 
    lighting gouraud
end



%% Get fiber volume coordinates and nucleus CoG coordinates
indices     = find(surf_mask_stack);
% Flip x and y dimensions to image notation
[y, x, z]   = ind2sub(size(surf_mask_stack), indices);

x = x.*voxel_size(1).*subsampling_factor;
y = y.*voxel_size(2).*subsampling_factor;
z = z.*voxel_size(3).*subsampling_factor;

voxel_loc_fiber = [x y z];



% Prepare results matrices used to store results for each CoG
% All distance measures are in um, MNDS in um^2, MNDV in um^3
% CoGs_inside_results:  [  Increment           |       Slice number  |   Slice center |
%                        Distance to fiber tip |               MNDS  |           MNDV |
%                           Inside CoG NN dist]
%
% CoGs_outside_results: [  Increment           |     Slice number(*) |       Slice center |
%                        Distance to fiber tip | Outside CoG NN dist | Inside CoG NN dist | Fiber surface vertex NN dist]
%
% (*) For negative values, the slice number is negative, increasing from 0 
%
CoGs_inside_results     = zeros(size(CoGs_inside, 1), 7);
CoGs_outside_results    = zeros(size(CoGs_outside, 1), 7);


if show_interm_figures
    hold on;
%     plot3(voxel_loc_fiber(1:10:end, 1), ...
%           voxel_loc_fiber(1:10:end, 2), ...
%           voxel_loc_fiber(1:10:end, 3), 'r.');
    plot3(CoGs_inside(:, 1), CoGs_inside(:, 2), CoGs_inside(:, 3), 'k*');
    plot3(CoGs_inside(:, 1), CoGs_inside(:, 2), CoGs_inside(:, 3), 'bo');
    plot3(CoGs_outside(:, 1), CoGs_outside(:, 2), CoGs_outside(:, 3), 'k*');
    plot3(CoGs_outside(:, 1), CoGs_outside(:, 2), CoGs_outside(:, 3), 'ro');
end

% Legacy code...
smoothed_nfv = fv;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use PCA to find the principal components and rotate the sample with the 
% first principal component along the z-axis and the second along the
% x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[v_pca, mean_pc, u, v] = rotatePointCloudAlongZ(smoothed_nfv.vertices, 'x');

% Rotational angle between the projected u on the xy plane and the x-axis
% (in Imaris)
alpha   = atan2(u(2), u(1)) * 180 / pi; 
%beta    = atan2(sqrt(u(1)^2 + u(2)^2), u(3)) * 180 / pi;

% Generate rotated mesh
fv_aligned          = smoothed_nfv;
fv_aligned.vertices = v_pca;

% Apply the rotation to the CoGs as well ...
CoGs_inside_pca  = rotatePointCloudAlongZ(CoGs_inside, 'x', mean_pc, u, v);
CoGs_outside_pca = rotatePointCloudAlongZ(CoGs_outside, 'x', mean_pc, u, v);

% ... and apply the rotation to the voxel volume of the fiber as well
voxel_loc_fiber_pca = rotatePointCloudAlongZ(voxel_loc_fiber, 'x', mean_pc, u, v);

% ... and finally move fiber surface, CoGs and the volume to the yz-plane,
% i.e. the tip of the fibre is located at x = 0
min_fv_aligned_vertices   = min(fv_aligned.vertices(:,3));
fv_aligned.vertices(:,3)  = fv_aligned.vertices(:,3) - min_fv_aligned_vertices;
CoGs_inside_pca(:,3)      = CoGs_inside_pca(:,3) - min_fv_aligned_vertices;
CoGs_outside_pca(:,3)     = CoGs_outside_pca(:,3) - min_fv_aligned_vertices;
voxel_loc_fiber_pca(:,3)  = voxel_loc_fiber_pca(:,3) - min_fv_aligned_vertices;

% Sometimes the data is 180 degrees wrongly rotated, so ask user to check
% and redo if necessary
surface_properly_aligned = 0;
use_dialog               = 1;

while ~surface_properly_aligned
    % Plot rotated mesh
    plot_surface_mesh(fv_aligned, 0.6);
    hold on;
    plot3(CoGs_inside_pca(:, 3), CoGs_inside_pca(:, 1), CoGs_inside_pca(:, 2), 'k*');
    plot3(CoGs_inside_pca(:, 3), CoGs_inside_pca(:, 1), CoGs_inside_pca(:, 2), 'bo');
    plot3(CoGs_outside_pca(:, 3), CoGs_outside_pca(:, 1), CoGs_outside_pca(:, 2), 'k*');
    plot3(CoGs_outside_pca(:, 3), CoGs_outside_pca(:, 1), CoGs_outside_pca(:, 2), 'ro');
    xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [\mum]');
    set(gcf, 'Position', [145         118        1537         795]);

    if use_dialog
        answer_align = questdlg('Data properly aligned?', 'Alignment', 'Yes','No', 'Yes');
        
        % Handle response
        if strcmp(answer_align, 'No')
            close gcf;        
            fv_aligned          = transform_surface(fv_aligned, get_rotation_matrix([0 180 0]));
            CoGs_inside_pca     = transform_pointset(CoGs_inside_pca, get_rotation_matrix([0 180 0]));
            CoGs_outside_pca    = transform_pointset(CoGs_outside_pca, get_rotation_matrix([0 180 0]));
            voxel_loc_fiber_pca = transform_pointset(voxel_loc_fiber_pca, get_rotation_matrix([0 180 0]));
        
            % ... and move the data to the yz-plane
            min_fv_aligned_vertices   = min(fv_aligned.vertices(:,3));
            fv_aligned.vertices(:,3)  = fv_aligned.vertices(:,3) - min_fv_aligned_vertices;
            CoGs_inside_pca(:,3)      = CoGs_inside_pca(:,3) - min_fv_aligned_vertices;
            CoGs_outside_pca(:,3)     = CoGs_outside_pca(:,3) - min_fv_aligned_vertices;
            voxel_loc_fiber_pca(:,3)  = voxel_loc_fiber_pca(:,3) - min_fv_aligned_vertices;
    
            use_dialog = 0;
        else
            surface_properly_aligned = 1;
        end
    else
        surface_properly_aligned = 1;
    end
end

h_surface_with_nuclei = gcf;

if debugging
    figure('Name', 'Surface and spots aligned');
    pause(.1);
else
    disp('Surface and spots aligned ...');
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate NN distances from nuclei to other nuclei and fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate distances between nuclei...');

dist_inside_nuclei              = [];
dist_outside_nuclei             = [];
dist_outside_inside_nuclei      = [];
nucleus_inside_indices          = 1:size(CoGs_inside_pca, 1);
nucleus_outside_indices         = 1:size(CoGs_outside_pca, 1);

% Inside NN for each inside nucleus
for nucleus_cnt = 1:size(CoGs_inside_pca, 1)
    tmp = repmat(CoGs_inside_pca(nucleus_cnt, 1:3), size(CoGs_inside_pca, 1) - 1, 1);
    
    % Exclude the current CoG from the calculation
    dist_tmp = ((tmp(:,1) - CoGs_inside_pca(find(nucleus_inside_indices ~= nucleus_cnt),1)).^2 + ...
                (tmp(:,2) - CoGs_inside_pca(find(nucleus_inside_indices ~= nucleus_cnt),2)).^2 + ...
                (tmp(:,3) - CoGs_inside_pca(find(nucleus_inside_indices ~= nucleus_cnt),3)).^2).^0.5;

    dist_inside_nuclei = [dist_inside_nuclei dist_tmp];
end

nucleus_inside_NNs       = min(dist_inside_nuclei);
mean_nucleus_inside_NNs  = mean(nucleus_inside_NNs);
std_nucleus_inside_NNs   = std(nucleus_inside_NNs);

% Fill CoG inside results matrix
CoGs_inside_results(:,1) = (1:size(CoGs_inside, 1))';
CoGs_inside_results(:,7) = nucleus_inside_NNs';



% Outside NN for each outside nucleus
for nucleus_cnt = 1:size(CoGs_outside_pca, 1)
    tmp = repmat(CoGs_outside_pca(nucleus_cnt, 1:3), size(CoGs_outside_pca, 1) - 1, 1);
    
    % Exclude the current CoG from the calculation
    dist_tmp = ((tmp(:,1) - CoGs_outside_pca(find(nucleus_outside_indices ~= nucleus_cnt),1)).^2 + ...
                (tmp(:,2) - CoGs_outside_pca(find(nucleus_outside_indices ~= nucleus_cnt),2)).^2 + ...
                (tmp(:,3) - CoGs_outside_pca(find(nucleus_outside_indices ~= nucleus_cnt),3)).^2).^0.5;

    dist_outside_nuclei = [dist_outside_nuclei dist_tmp];
end

nucleus_outside_NNs       = min(dist_outside_nuclei);
mean_nucleus_outside_NNs  = mean(nucleus_outside_NNs);
std_nucleus_outside_NNs   = std(nucleus_outside_NNs);

% Fill CoG outside results matrix
CoGs_outside_results(:,1) = 1:size(CoGs_outside, 1);
CoGs_outside_results(:,5) = nucleus_outside_NNs;



% Inside NN for each outside nucleus
for nucleus_cnt = 1:size(CoGs_outside_pca, 1)
    tmp = repmat(CoGs_outside_pca(nucleus_cnt, 1:3), size(CoGs_inside_pca, 1), 1);
    
    % Exclude the current CoG from the calculation
    dist_tmp = ((tmp(:,1) - CoGs_inside_pca(:,1)).^2 + ...
                (tmp(:,2) - CoGs_inside_pca(:,2)).^2 + ...
                (tmp(:,3) - CoGs_inside_pca(:,3)).^2).^0.5;

    dist_outside_inside_nuclei = [dist_outside_inside_nuclei dist_tmp];
end

nucleus_outside_inside_NNs      = min(dist_outside_inside_nuclei);
mean_nucleus_outside_inside_NNs = mean(nucleus_outside_inside_NNs);
std_nucleus_outside_inside_NNs  = std(nucleus_outside_inside_NNs);

% Fill CoG outside results matrix
CoGs_outside_results(:,6) = nucleus_outside_inside_NNs;



%% Calculate NN distances from outside nuclei to fibre surface vertices
disp('Calculate NN distances between outside nuclei and surface vertices...');

dist_outside_nuclei_surface     = [];

% Surface vertex NN for each outside nucleus
for nucleus_cnt = 1:size(CoGs_outside_pca, 1)
    tmp = repmat(CoGs_outside_pca(nucleus_cnt, 1:3), size(fv_aligned.vertices, 1), 1);
    
    % Exclude the current CoG from the calculation
    dist_tmp = ((tmp(:,1) - fv_aligned.vertices(:,1)).^2 + ...
                (tmp(:,2) - fv_aligned.vertices(:,2)).^2 + ...
                (tmp(:,3) - fv_aligned.vertices(:,3)).^2).^0.5;

    dist_outside_nuclei_surface = [dist_outside_nuclei_surface dist_tmp];
end

nucleus_outside_surface_NNs      = min(dist_outside_nuclei_surface);
mean_nucleus_outside_surface_NNs = mean(nucleus_outside_surface_NNs);
std_nucleus_outside_surface_NNs  = std(nucleus_outside_surface_NNs);

% Fill CoG outside results matrix
CoGs_outside_results(:,7)        = nucleus_outside_surface_NNs;



%% Calculate NN distances from fiber surface vertices to inside CoGs. This 
% will later be necessary to calculate which surface vertices and faces
% are part of the myonucleus domain
disp('Calculate NN distances between surface vertices and inside nuclei...');
dist_to_surf = [];

for nucleus_cnt = 1:size(CoGs_inside_pca, 1)
    tmp = repmat(CoGs_inside_pca(nucleus_cnt, :), size(fv_aligned.vertices, 1), 1);
    
    distance_tmp = ((tmp(:,1) - fv_aligned.vertices(:,1)).^2 + ...
                    (tmp(:,2) - fv_aligned.vertices(:,2)).^2 + ...
                    (tmp(:,3) - fv_aligned.vertices(:,3)).^2).^0.5;

    dist_to_surf = [dist_to_surf distance_tmp];
end

% Find the nucleus with the shortest distance to each surface vertex
[~, idx]             = min(dist_to_surf');

fv_aligned_vertices  = fv_aligned.vertices;
fv_aligned_faces     = fv_aligned.faces;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the MNDS, thus the surface area belonging to each nucleus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate MND surface area...');

% Generate face area list
vectors_1st = fv_aligned_vertices(fv.faces(:,1), :);
vectors_2nd = fv_aligned_vertices(fv.faces(:,2), :);
vectors_3rd = fv_aligned_vertices(fv.faces(:,3), :);
        
% Vector cross product
tmp         = cross(vectors_3rd - vectors_1st, vectors_2nd - vectors_1st); 

% Half the vector norm is the triangle area
face_areas  = ((tmp(:,1).^2 + tmp(:,2).^2 + tmp(:,3).^2).^0.5)/2;

% Assign a face area to each face and determine an average face area for each
% vertex
face_areas                   = [face_areas face_areas face_areas]';
fv_faces                     = fv.faces';
average_face_area_per_vertex = accumarray(fv_faces(:), face_areas(:), [], @(x) mean(x,1)).*2;

% Determine the vertices belonging to each MNDS and calculate the surface
% that belongs to the patch
surface_per_label            = zeros(size(CoGs_inside_pca, 1), 1);

for nucleus_cnt = 1:size(CoGs_inside_pca, 1)
    cur_label_indices = find(idx == nucleus_cnt);   % Vertices closest to the current nucleus

    % Check if the current patch is at the end of the fiber (max. x-value in
    % Imaris). If so, exclude.
    if max(fv_orig(cur_label_indices, 1)) == size(surf_mask_stack, 2)
        surface_per_label(nucleus_cnt) = nan;
    else
        surface_per_label(nucleus_cnt) = sum(average_face_area_per_vertex(cur_label_indices));
    end        
end

% Calculate the mean and std for the patches, NaN 'resistent'
mean_surface_per_label = nanmean(surface_per_label);
std_surface_per_label  = nanstd(surface_per_label);

clear cur_label_indices_tmp;
clear fv_aligned_faces_tmp;
clear fv_aligned_faces_2;

% Fill CoG inside results matrix
CoGs_inside_results(:, 5) = surface_per_label;

% Prepare results plot. First generate face colors...
color_increment = floor(255 / size(CoGs_inside_pca, 1));
colorfield      = idx';
colorfield      = (colorfield * color_increment) - (color_increment - 1);      

edgecolor   = 'none';

% ... then plot
h_surface_opaque = figure;
colormap(jet);
trimesh(fv_aligned_faces, fv_aligned_vertices(:,3), fv_aligned_vertices(:,1), fv_aligned_vertices(:,2), ...
                colorfield, ...
                'FaceColor', 'flat', ...
                'EdgeColor', edgecolor, 'FaceAlpha', 1);

daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
set(gcf, 'Position', [145         118        1537         795]);


h_surface = figure;
colormap(jet);
trimesh(fv_aligned_faces, fv_aligned_vertices(:,3), fv_aligned_vertices(:,1), fv_aligned_vertices(:,2), ...
                colorfield, ...
                'FaceColor', 'flat', ...
                'EdgeColor', edgecolor, 'FaceAlpha', 0.4);

daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
set(gcf, 'Position', [145         118        1537         795]);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the MNDV, thus the volume belonging to each nucleus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate shortest distances from fiber voxels to nuclei...');
dist_to_vol     = [];

for nucleus_cnt = 1:size(CoGs_inside_pca, 1)
    tmp = repmat(CoGs_inside_pca(nucleus_cnt, 1:3), size(voxel_loc_fiber_pca, 1), 1);
    
    dist_tmp = ((tmp(:,1) - voxel_loc_fiber_pca(:,1)).^2 + ...
                (tmp(:,2) - voxel_loc_fiber_pca(:,2)).^2 + ...
                (tmp(:,3) - voxel_loc_fiber_pca(:,3)).^2).^0.5;

    dist_to_vol = [dist_to_vol dist_tmp];
end

% Find the nucleus with the shortest distance to each surface vertex
[~, idx] = min(dist_to_vol');
idx      = idx';

fiber_volume        = length(indices) * voxel_volume;
volumes_per_label   = [];

% Get fiber mask coordinates in x (Imaris)
indices                 = find(surf_mask_stack);
% Flip x and y dimensions to image notation
[y_tmp, x_tmp, z_tmp]   = ind2sub(size(surf_mask_stack), indices);

for nucleus_cnt = 1:size(CoGs_inside_pca, 1)
    % Check if the current volume is at the end of the fiber (max. x-value in
    % Imaris). If so, exclude.
    cur_label_volume = [x_tmp(idx == nucleus_cnt) y_tmp(idx == nucleus_cnt) z_tmp(idx == nucleus_cnt)];

    if max(cur_label_volume(:, 1)) == size(surf_mask_stack, 2)
        volumes_per_label = [volumes_per_label; nan];
    else
        volumes_per_label = [volumes_per_label; sum(idx == nucleus_cnt)];
    end
end

% Calculate the mean and std for the volumes, NaN 'resistent'
volumes_per_label       = volumes_per_label.*voxel_volume;
mean_volumes_per_label  = nanmean(volumes_per_label);
std_volumes_per_label   = nanstd(volumes_per_label);

% Fill CoG inside results matrix
CoGs_inside_results(:, 6) = volumes_per_label;


% scatter3 does not accept a blob size of 0, so use a very small value...
volumes_per_label_plotting                          = volumes_per_label;
volumes_per_label_plotting(volumes_per_label == 0)  = 0.00000001;


hold on;
c = colormap('jet');

for CoG_cnt = 1:size(CoGs_inside_pca, 1)
    scatter3(CoGs_inside_pca(CoG_cnt, 3), ...
             CoGs_inside_pca(CoG_cnt, 1), ...
             CoGs_inside_pca(CoG_cnt, 2), ...
             volumes_per_label_plotting(CoG_cnt)./5, ...
             c((CoG_cnt - 1) * color_increment + color_increment, :), 'filled');
    drawnow;
end

daspect([1 1 1]);
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [\mum]');




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subdivide the volume, determine the number of inside CoGs and the number 
%  of outside CoGs per slice. In addition, link the MNDV, MNDS as well as 
%  inside to inside nuclei, outside to outside nuclei, outside to inside 
%  nuclei and outside to surface nuclei to each slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Determine nucleus count per slice, slice nucleus density and ' ...
      'nucleus to nucleus distances...']);

% Get the last quarter of the fiber along the longitudinal axis ...
voxel_loc_fiber_pca_last_3rd = ...
    voxel_loc_fiber_pca(voxel_loc_fiber_pca(:, 3) > max(voxel_loc_fiber_pca(:,3))*0.75, :);

% ... project it to the yz-plane (Imaris reference) and get the maximum
% bounding box dimension of the shape
max_diameter_last_3rd = max(max(voxel_loc_fiber_pca_last_3rd(:,[1 2])) - ...
                        min(voxel_loc_fiber_pca_last_3rd(:,[1 2])));

% Now determine, how many slices have to be removed from analysis at the far 
% end of the fiber (bc then might not contain data in the entire slice)
cut_off_at_end        = sind(alpha) * max_diameter_last_3rd;
cut_off_at_end_slices = ceil(abs(cut_off_at_end / slice_thickness));

num_slices_inside  = ceil(max(fv_aligned.vertices(:,3))/slice_thickness) - ...
                     cut_off_at_end_slices;
slice_edges_inside = 0:slice_thickness:num_slices_inside * slice_thickness;


inside_nuclei_per_slice        = [];
outside_nuclei_per_slice       = [];

% In this cell array, each row corresponds to a slice. The columns contain 
% the following:
%
% 1st column                | 2nd column                | 3rd column                | 4th column
% Indices to the nuclei     | Inside NN distance        | Surface area              | Volume 
inside_nuclei_labels_per_slice  = cell(num_slices_inside, 4);

volume_per_slice        = [];

CoG_inside_z_values     = CoGs_inside_pca(:,3);
CoG_outside_z_values    = CoGs_outside_pca(:,3);
voxel_loc_z_values      = voxel_loc_fiber_pca(:,3);


% Determine slice inside CoG count and slice volumes
for cur_slice = 1:num_slices_inside
    % Find CoGs in current slice
    cond_1 = CoG_inside_z_values >= slice_edges_inside(cur_slice);
    cond_2 = CoG_inside_z_values < slice_edges_inside(cur_slice + 1);

    inside_nuclei_per_slice          = [inside_nuclei_per_slice; sum(cond_1 & cond_2)];
    volume_per_slice                 = [volume_per_slice; ...
                                        sum(voxel_loc_z_values >= slice_edges_inside(cur_slice) & ...
                                        voxel_loc_z_values < slice_edges_inside(cur_slice + 1)).*voxel_volume];    

    tmp_indices                                  = find(cond_1 & cond_2);
    inside_nuclei_labels_per_slice{cur_slice, 1} = tmp_indices;
    inside_nuclei_labels_per_slice{cur_slice, 2} = nucleus_inside_NNs(tmp_indices);

    % Fill CoG inside results matrix
    CoGs_inside_results(tmp_indices, 2)          = cur_slice;
    CoGs_inside_results(tmp_indices, 3)          = (slice_edges_inside(cur_slice) + slice_edges_inside(cur_slice + 1))/2;

    % Also add surface and volume information for each inside nucleus
    if ~isempty(tmp_indices)
        inside_nuclei_labels_per_slice{cur_slice, 3} = surface_per_label(tmp_indices);
        inside_nuclei_labels_per_slice{cur_slice, 4} = volumes_per_label(tmp_indices);
    end
end

% Fill CoG outside results matrix
% These are in slices that touch the border of the image domain
CoGs_inside_results(CoGs_inside_results(:, 2) == 0, 2) = nan;
CoGs_inside_results(CoGs_inside_results(:, 2) == 0, 3) = nan;
% Add distance to tip
CoGs_inside_results(:, 4)                              = CoG_inside_z_values;



% Determine slice outside CoG count
if min(CoG_outside_z_values) < 0
    num_neg_slices = abs(floor(min(CoG_outside_z_values)/slice_thickness));
else
    num_neg_slices = 0;
end

num_pos_slices      = ceil(max(fv_aligned.vertices(:,3))/slice_thickness) - ...
                      cut_off_at_end_slices;
slice_edges_outside = -(num_neg_slices * slice_thickness):slice_thickness:num_pos_slices * slice_thickness;
zero_edge_slice     = find(slice_edges_outside == 0);

% In this cell array, each row corresponds to a slice. The columns contain 
% the following:
%
% 1st column                | 2nd column                | 3rd column                | 4th column
% Indices to the nuclei     | Outside NN distance       | Outside to inside NN dist | Outside to NN surface vertex distance
outside_nuclei_labels_per_slice = cell(length(slice_edges_outside) - 1, 4);

for cur_slice = 1:length(slice_edges_outside) - 1
    % Find CoGs in current slice
    cond_1 = CoG_outside_z_values >= slice_edges_outside(cur_slice);
    cond_2 = CoG_outside_z_values < slice_edges_outside(cur_slice + 1);

    outside_nuclei_per_slice                      = [outside_nuclei_per_slice; sum(cond_1 & cond_2)];

    tmp_indices                                   = find(cond_1 & cond_2);
    outside_nuclei_labels_per_slice{cur_slice, 1} = tmp_indices;
    outside_nuclei_labels_per_slice{cur_slice, 2} = nucleus_outside_NNs(tmp_indices);

    % Fill CoG outside results matrix
    % For each CoG with a positive distance, the slices are counted 1, 2
    % ... n, for negative distances, the slices are counted -1, -2 ... -m
    if cur_slice < zero_edge_slice
        CoGs_outside_results(tmp_indices, 2) = -(zero_edge_slice - cur_slice);
    else
        CoGs_outside_results(tmp_indices, 2) = cur_slice - zero_edge_slice + 1;
    end

    CoGs_outside_results(tmp_indices, 3) = (slice_edges_outside(cur_slice) + slice_edges_outside(cur_slice + 1))/2;

    % Also add outside to inside nucleus NN distance as well as outside to
    % fiber surface NN vertex for each outside nucleus
    if ~isempty(tmp_indices)
        outside_nuclei_labels_per_slice{cur_slice, 3} = nucleus_outside_inside_NNs(tmp_indices);
        outside_nuclei_labels_per_slice{cur_slice, 4} = nucleus_outside_surface_NNs(tmp_indices);
    end    
end

% Fill CoG outside results matrix
% These are in slices that touch the border of the image domain
CoGs_outside_results(CoGs_outside_results(:, 2) == 0, 2) = nan;
CoGs_outside_results(CoGs_outside_results(:, 2) == 0, 3) = nan;
% Add distance to tip
CoGs_outside_results(:, 4)                               = CoG_outside_z_values;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot and save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if results folder for the current results exists. If so, ask the
% user whether to overwrite results or not. If a results folder cannot be
% created by one or the other reason, the results can still be studied but
% are not saved automatically.
save_results = 0;

if ~exist(fullfile(results_saving_folder, data_ID{1}), 'dir')
    [success, ~, ~] = mkdir(fullfile(results_saving_folder, data_ID{1}));

    if ~success
        answer_no_save = questdlg('Results will not be saved, continue?', 'Save results', 'Yes','No', 'Yes');

        if strcmp(answer_no_save, 'No')
            return;
        end
    else
        save_results = 1;
    end
else
    answer_save = questdlg('Results exist, overwrite?', 'Save results', 'Yes','No', 'Yes');
    
    % Handle response
    if strcmp(answer_save, 'Yes')
        [success, ~, ~] = rmdir(fullfile(results_saving_folder, data_ID{1}), 's');
        [success, ~, ~] = mkdir(fullfile(results_saving_folder, data_ID{1}));

        if ~success
            answer_no_save = questdlg('Results will not be saved, continue?', 'Save results', 'Yes','No', 'Yes');

            if strcmp(answer_no_save, 'No')
                return;
            end
        else
            save_results = 1;
        end
    end
end

if save_results
    savefig(h_surface, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_MNDS_MNDV.fig'));
    saveas(h_surface, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_MNDS_MNDV.png'));
    savefig(h_surface_opaque, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_MNDS_MNDV_opaque.fig'));
    saveas(h_surface_opaque, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_MNDS_MNDV_opaque.png'));
    savefig(h_surface_with_nuclei, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_nuclei_in_out.fig'));
    saveas(h_surface_with_nuclei, fullfile(results_saving_folder, data_ID{1}, 'Fiber_surf_nuclei_in_out.png'));
end



% Plot number of inside nuclei per slice
figure;
subplot(1,2,1);
plot(slice_edges_inside(2:end) - slice_thickness/2, inside_nuclei_per_slice, 'b*-');
hold on;
plot(slice_edges_outside(2:end) - slice_thickness/2, outside_nuclei_per_slice, 'r*-');
title(['Nuclei per slice (thickness = ' num2str(slice_thickness) '\mum)']);
xlabel('Slice center location [\mum]');
ylabel('# of nuclei');
grid on;
legend('Inside nuclei', 'Outside nuclei');
x_lim = get(gca, 'XLim');
x_lim(1) = x_lim(1) - 10;
x_lim(2) = x_lim(2) + 10;
set(gca, 'XLim', x_lim);

% Plot number of inside nuclei density per slice
subplot(1,2,2);
plot(slice_edges_inside(2:end) - slice_thickness/2, inside_nuclei_per_slice./volume_per_slice, 'b*-');
title('Nuclei density per slice');
xlabel('Slice center location [\mum]');
ylabel('# of nuclei per volume [\mum^3]');
grid on;
set(gca, 'XLim', x_lim);

set(gcf, 'Position', [250         557        1026         371]);

if save_results
    savefig(gcf, fullfile(results_saving_folder, data_ID{1}, 'number_density_nuclei_per_slice.fig'));
    saveas(gcf, fullfile(results_saving_folder, data_ID{1}, 'number_density_nuclei_per_slice.png'));
end



% Plot surface area and volume per slice
figure;
subplot(1,2,1);
hold on;
for cur_slice = 1:num_slices_inside
    if length(inside_nuclei_labels_per_slice{cur_slice, 3}) > 0
        for cur_nucleus = 1:length(inside_nuclei_labels_per_slice{cur_slice, 3})
            plot(slice_edges_inside(cur_slice + 1) - slice_thickness/2, inside_nuclei_labels_per_slice{cur_slice, 3}(cur_nucleus), 'k*');
        end
    end
end
title(['MNDS (\mum^2), ' num2str(mean_surface_per_label) ' +/- ' num2str(std_surface_per_label)]);
xlabel('Slice center location [\mum]');
ylabel('Myonuclei domain surface [\mum^2]');
set(gca, 'Xlim', x_lim);
grid on;


subplot(1,2,2);
hold on;
for cur_slice = 1:num_slices_inside
    if length(inside_nuclei_labels_per_slice{cur_slice, 3}) > 0
        for cur_nucleus = 1:length(inside_nuclei_labels_per_slice{cur_slice, 3})
            plot(slice_edges_inside(cur_slice + 1) - slice_thickness/2, inside_nuclei_labels_per_slice{cur_slice, 4}(cur_nucleus), 'c*');
        end
    end
end
title(['MNDV (\mum^3), ' num2str(mean_volumes_per_label) ' +/- ' num2str(std_volumes_per_label)]);
xlabel('Slice center location [\mum]');
ylabel('Myonuclei domain volume [\mum^3]');
set(gca, 'Xlim', x_lim);
grid on;

set(gcf, 'Position', [250         557        1026         371]);

if save_results
    savefig(gcf, fullfile(results_saving_folder, data_ID{1}, 'MNDS_MNDV_per_slice.fig'));
    saveas(gcf, fullfile(results_saving_folder, data_ID{1}, 'MNDS_MNDV_per_slice.png'));
end



% Plot inside and outside nucleus NN distance per slice
figure;
subplot(1,3,1);
hold on;
for cur_slice = 1:num_slices_inside
    if length(inside_nuclei_labels_per_slice{cur_slice, 2}) > 0
        for cur_nucleus = 1:length(inside_nuclei_labels_per_slice{cur_slice, 2})
            plot(slice_edges_inside(cur_slice + 1) - slice_thickness/2, inside_nuclei_labels_per_slice{cur_slice, 2}(cur_nucleus), 'b*');
        end
    end
end
for cur_slice = 1:length(slice_edges_outside) - 1
    if length(outside_nuclei_labels_per_slice{cur_slice, 2}) > 0
        for cur_nucleus = 1:length(outside_nuclei_labels_per_slice{cur_slice, 2})
            plot(slice_edges_outside(cur_slice + 1) - slice_thickness/2, outside_nuclei_labels_per_slice{cur_slice, 2}(cur_nucleus), 'r*');
        end
    end
end
title(['Inside (blue, \mum), ' num2str(mean_nucleus_inside_NNs) ' +/- ' num2str(std_nucleus_inside_NNs), ...
       ', Outside (red, \mum), ' num2str(mean_nucleus_outside_NNs) ' +/- ' num2str(std_nucleus_outside_NNs)]);
xlabel('Slice center location [\mum]');
ylabel('Distance to nearest nucleus [\mum]');
set(gca, 'Xlim', x_lim);
grid on;


% Plot outside to inside nucleus NN distance per slice
subplot(1,3,2);
hold on;
for cur_slice = 1:length(slice_edges_outside) - 1
    if length(outside_nuclei_labels_per_slice{cur_slice, 3}) > 0
        for cur_nucleus = 1:length(outside_nuclei_labels_per_slice{cur_slice, 3})
            plot(slice_edges_outside(cur_slice + 1) - slice_thickness/2, ...
                 outside_nuclei_labels_per_slice{cur_slice, 3}(cur_nucleus), ...
                 'Color', [1 0 1], 'Marker', '*');
        end
    end
end
title(['Outside to inside (\mum), ' num2str(mean_nucleus_outside_inside_NNs) ' +/- ' num2str(std_nucleus_outside_inside_NNs)]);
xlabel('Slice center location [\mum]');
ylabel('Distance to nearest inside nucleus [\mum]');
set(gca, 'Xlim', x_lim);
grid on;


% Plot outside to surface vertex NN distance per slice
subplot(1,3,3);
hold on;
for cur_slice = 1:length(slice_edges_outside) - 1
    if length(outside_nuclei_labels_per_slice{cur_slice, 4}) > 0
        for cur_nucleus = 1:length(outside_nuclei_labels_per_slice{cur_slice, 4})
            plot(slice_edges_outside(cur_slice + 1) - slice_thickness/2, ...
                 outside_nuclei_labels_per_slice{cur_slice, 4}(cur_nucleus), ...
                 'Color', [1 0.75 0], 'Marker', '*');
        end
    end
end
title(['Outside to surface (\mum), ' num2str(mean_nucleus_outside_surface_NNs) ' +/- ' num2str(std_nucleus_outside_surface_NNs)]);
xlabel('Slice center location [\mum]');
ylabel('Distance to nearest fiber surface vertex [\mum]');
set(gca, 'Xlim', x_lim);
grid on;

set(gcf, 'Position', [139         159        1578         371]);

if save_results
    savefig(gcf, fullfile(results_saving_folder, data_ID{1}, 'inter_nuclei_and_to_surface_distance_per_slice.fig'));
    saveas(gcf, fullfile(results_saving_folder, data_ID{1}, 'inter_nuclei_and_to_surface_distance_per_slice.png'));
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a labeled fiber volume for visualization and save it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeled_fiber = zeros(size(surf_mask_stack, 1), size(surf_mask_stack, 2), size(surf_mask_stack, 3), 'uint8');

indices     = find(surf_mask_stack);
[x, y, z]   = ind2sub(size(surf_mask_stack), indices);

[~, idx] = min(dist_to_vol');
idx      = idx';

for cur_vox = 1:size(idx, 1)
    labeled_fiber(x(cur_vox), y(cur_vox), z(cur_vox)) = idx(cur_vox);
end

orthosliceViewer(labeled_fiber.*floor(255/max(labeled_fiber(:))));
%dipshow(labeled_fiber.*floor(255/max(labeled_fiber(:))));

if save_results
    [success, ~, ~] = rmdir(fullfile(results_saving_folder, data_ID{1}, 'labeled_volume'), 's');
    [success, ~, ~] = mkdir(fullfile(results_saving_folder, data_ID{1}, 'labeled_volume'));

    if success
        tmp = labeled_fiber.*floor(255/max(labeled_fiber(:)));

        for image_cnt = 1:size(tmp, 3)
            imwrite(tmp(:,:,image_cnt), fullfile(results_saving_folder, data_ID{1}, ...
                   'labeled_volume', ['Slice_' num2str(image_cnt) '.png']));    
        end
    end    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export results to Excel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_results
    if exist(excel_file_path, 'file')
        % Get CoG inside data, remove 'missing' cells and determine starting
        % row for current block. If current data_ID already exists, overwrite
        % that data block
        xls_inside = readcell(excel_file_path, 'Sheet', 'Inside nuclei data');
        [xls_inside, starting_row_CoGs_inside_results, ending_row_CoGs_inside_results] = ...
                        remove_missing_cells_determine_starting_row(xls_inside, data_ID);
    
        % Get CoG outside data and ...
        xls_outside = readcell(excel_file_path, 'Sheet', 'Outside nuclei data');
        [xls_outside, starting_row_CoGs_outside_results, ending_row_CoGs_outside_results] = ...
                        remove_missing_cells_determine_starting_row(xls_outside, data_ID);
    
        % Get inside slice data and ...
        xls_inside_slice = readcell(excel_file_path, 'Sheet', 'Inside slice data');
        [xls_inside_slice, starting_row_inside_slice_results, ~] = ...
                        remove_missing_cells_determine_starting_row(xls_inside_slice, data_ID);
    
        % Get outside slice data and ...
        xls_outside_slice = readcell(excel_file_path, 'Sheet', 'Outside slice data');
        [xls_outside_slice, starting_row_outside_slice_results, ~] = ...
                        remove_missing_cells_determine_starting_row(xls_outside_slice, data_ID);
    else
        slices_edges_positive = num2cell(0:slice_thickness:1500);   % in um
        slices_edges_negative = num2cell(0:-slice_thickness:-1500); % in um
    
        xls_slices                                      = cell(6, size(slices_edges_positive, 2) + 1);
    
        xls_slices(1,1)                                 = {'Slice edges positive direction'};
        xls_slices(2,1)                                 = {'Slice count'};
        xls_slices(3,1)                                 = {'Edges'};
        xls_slices(2,3:size(slices_edges_positive, 2) + 1)  = num2cell(1:(size(slices_edges_positive, 2)-1));
        xls_slices(3,2:size(slices_edges_positive, 2) + 1)  = slices_edges_positive;
    
        xls_slices(4,1)                                 = {'Slice edges negative direction'};
        xls_slices(5,1)                                 = {'Slice count'};
        xls_slices(6,1)                                 = {'Edges'};    
        xls_slices(5,3:size(slices_edges_negative, 2) + 1)  = num2cell(-1:-1:-(size(slices_edges_negative, 2)-1));
        xls_slices(6,2:size(slices_edges_negative, 2) + 1)  = slices_edges_negative;
        
        writecell(xls_slices, excel_file_path, 'Sheet', 'Slices');
    
    
        xls_inside_slice        = cell(5, size(slice_edges_inside, 2) + 2);
        xls_inside_slice(1,1:3) = [{'Data IDs'} {'Label'} {'Slice data'}];
    
        xls_outside_slice        = cell(4, size(slice_edges_outside, 2) + 2);
        xls_outside_slice(1,1:3) = [{'Data IDs'} {'Label'} {'Slice data'}];
    
        xls_inside      = cell(size(CoGs_inside_results, 1) + 1, 8);
        xls_inside(1,:) = [{'Data IDs'} {'Count'} {'Slice number'} {'Slice center'} ...
                           {'Dist to fiber tip [um]'} {'MNDS [um^2]'} {'MNDV [um^3]'} ...
                           {'Inside-Inside NN [um]'}];
    
        xls_outside      = cell(size(CoGs_outside_results, 1) + 1, 8);
        xls_outside(1,:) = [{'Data IDs'} {'Count'} {'Slice number'} {'Slice center'} ...
                            {'Dist to fiber tip [um]'} {'Outside-Outside NN [um]'} ...
                            {'Outside-Inside NN [um]'} {'Outside-Fiber vertex NN [um]'}];
    
        starting_row_CoGs_inside_results   = 2;
        starting_row_CoGs_outside_results  = 2;
        starting_row_inside_slice_results  = 2;
        starting_row_outside_slice_results = 2;
    
        ending_row_CoGs_inside_results  = 2 + size(CoGs_inside_results, 1) - 1;
        ending_row_CoGs_outside_results = 2 + size(CoGs_outside_results, 1) - 1;
    end
    
    error_occurred = 0;
    
    % In case an Excel sheets exist and the current dataID was used before, the
    % program checks whether the size of the data block is the same as the
    % current data block to be saved. If yes, the data is overwritten. If not,
    % the dataID was used for another dataset before and in that case an error
    % is rendered to make sure the results file is not corrupted.
    if (ending_row_CoGs_inside_results - starting_row_CoGs_inside_results + 1) ~= ...
       (size(CoGs_inside_results, 1)) && ~isnan(ending_row_CoGs_inside_results) 
        f = errordlg(['Already existing data block in the Excel sheet with the current ID differs in size from the current result. ' ...
               'Could it be that this ID was used before for another dataset?'], ...
               'Possibly reuse of data ID ...' );
        error_occurred = 1;
    end
    
    if (ending_row_CoGs_outside_results - starting_row_CoGs_outside_results + 1) ~= ...
            (size(CoGs_outside_results, 1)) && ~isnan(ending_row_CoGs_outside_results) && ~error_occurred
        f = errordlg(['Already existing data block in the Excel sheet with the current ID differs in size from the current result. ' ...
               'Could it be that this ID was used before for another dataset?'], ...
               'Possibly reuse of data ID ...' );
        error_occurred = 1;    
    end
    
    if ~error_occurred
        cur_row                     = starting_row_CoGs_inside_results;
        xls_inside(cur_row:cur_row + size(CoGs_inside_results, 1) - 1, 1)     = data_ID;
        xls_inside(cur_row:cur_row + size(CoGs_inside_results, 1) - 1, 2:end) = num2cell(CoGs_inside_results);
    
        cur_row                     = starting_row_CoGs_outside_results;
        xls_outside(cur_row:cur_row + size(CoGs_outside_results, 1) - 1, 1)     = data_ID;
        xls_outside(cur_row:cur_row + size(CoGs_outside_results, 1) - 1, 2:end) = num2cell(CoGs_outside_results);
        
        writecell(xls_inside, excel_file_path, 'Sheet', 'Inside nuclei data');    
        writecell(xls_outside, excel_file_path, 'Sheet', 'Outside nuclei data');    
    
        cur_row                                   = starting_row_inside_slice_results;
        xls_inside_slice(cur_row:cur_row + 5, 1)  = data_ID;
        xls_inside_slice(cur_row, 2)              = {'Slice count'};
        xls_inside_slice(cur_row + 1, 2)          = {'Edges [um]'};
        xls_inside_slice(cur_row + 2, 2)          = {'Slice center [um]'};
        xls_inside_slice(cur_row + 3, 2)          = {'Nuclei per slice'};
        xls_inside_slice(cur_row + 4, 2)          = {'Slice volume [um^3]'};
        xls_inside_slice(cur_row + 5, 2)          = {'Slice nuclei density (nuclei per volume [um^3])'};
            
        xls_inside_slice(cur_row, 4:(4 + size(slice_edges_inside, 2) - 2))       = num2cell(1:(size(slice_edges_inside, 2) - 1));
        xls_inside_slice(cur_row + 1, 3:(3 + size(slice_edges_inside, 2) - 1))   = num2cell(slice_edges_inside);
        xls_inside_slice(cur_row + 2, 4:(4 + size(slice_edges_inside, 2) - 2))   = num2cell((slice_edges_inside(1:end-1) + slice_edges_inside(2:end))/2);
        xls_inside_slice(cur_row + 3, 4:(4 + size(slice_edges_inside, 2) - 2))   = num2cell(inside_nuclei_per_slice');
        xls_inside_slice(cur_row + 4, 4:(4 + size(slice_edges_inside, 2) - 2))   = num2cell(volume_per_slice');
        xls_inside_slice(cur_row + 5, 4:(4 + size(slice_edges_inside, 2) - 2))   = num2cell((inside_nuclei_per_slice./volume_per_slice)');
        
        
        cur_row                                     = starting_row_outside_slice_results;
        xls_outside_slice(cur_row:cur_row + 5, 1)   = data_ID;
        xls_outside_slice(cur_row, 2)               = {'Slice count'};
        xls_outside_slice(cur_row + 1, 2)           = {'Edges [um]'};
        xls_outside_slice(cur_row + 2, 2)           = {'Slice center [um]'};
        xls_outside_slice(cur_row + 3, 2)           = {'Nuclei per slice'};
        
        xls_outside_slice(cur_row, 4:(4 + size(slice_edges_outside, 2) - 2))     = num2cell(1:(size(slice_edges_outside, 2) - 1));
        xls_outside_slice(cur_row + 1, 3:(3 + size(slice_edges_outside, 2) - 1)) = num2cell(slice_edges_outside);
        xls_outside_slice(cur_row + 2, 4:(4 + size(slice_edges_outside, 2) - 2)) = num2cell((slice_edges_outside(1:end-1) + slice_edges_outside(2:end))/2);
        xls_outside_slice(cur_row + 3, 4:(4 + size(slice_edges_outside, 2) - 2)) = num2cell(outside_nuclei_per_slice');
        
        writecell(xls_inside_slice, excel_file_path, 'Sheet', 'Inside slice data');
        writecell(xls_outside_slice, excel_file_path, 'Sheet', 'Outside slice data');
    end
end



%% Ask to close all figures
answer_close_all = questdlg('Close all figures?', 'Close figures', 'Yes','No', 'Yes');

if strcmp(answer_close_all, 'Yes')
    close all;
end

end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_surface_mesh(fv, varargin)
    if nargin > 1
        face_alpha = varargin{1};
    else
        face_alpha = 1;
    end

    figure;
    colorfield_tmp = ones(length(fv.vertices(:,1)), 1);
    facecolor_tmp = [1 0.8 0.4];
    edgecolor_tmp = 'none';
    trimesh(fv.faces, fv.vertices(:,3), fv.vertices(:,1), ...
                      fv.vertices(:,2), ...
                      colorfield_tmp, 'FaceColor', facecolor_tmp, ...
                     'EdgeColor', edgecolor_tmp, 'FaceAlpha', face_alpha);
    daspect([1 1 1])
    view(3); axis tight
    camlight 
    lighting gouraud
end


function [output, starting_row, ending_row] = remove_missing_cells_determine_starting_row(input, data_ID)
    output = input;

    % Remove 'missing' cells ...
    for cell_cnt = 1:size(input(:), 1)
        if ismissing(input{cell_cnt})
            output{cell_cnt} = [];
        end
    end

    % ... and determine whether the current data_ID already exists. If so,
    % overwrite that data block
    ids          = output(2:end, 1);
    starting_row = size(output, 1) + 1;

    for ids_cnt = 1:length(ids)
        if strcmp(ids{ids_cnt}, data_ID{1})
            starting_row = ids_cnt + 1;
            break;
        end
    end

    % Also check the block size in the file. This is not really necessary,
    % but the user might forget to change the ID for the current dataset
    % and then that datablock (and possibly following ones) get
    % corrupted... so let's be sure
    ending_row = nan;
    if starting_row ~= (size(input, 1) + 1) 
        for ids_cnt = starting_row + 1:length(ids)
            if ~strcmp(ids{ids_cnt}, data_ID{1})
                ending_row = ids_cnt;
                break;
            end
        end        
    end
end


%%
function data_homo = cart2homo(data_eucl)
    % Converts a dataset in euclidean space (nx3) to homogeneous space (nx4)
    data_homo = [data_eucl ones(length(data_eucl(:,1)), 1)];
end


function data_eucl = homo2cart(data_homo)
    % Converts a dataset in homogeneous space to euclidean space
    data_eucl = data_homo(:, 1:3);
end


function rot_matrix  = get_rotation_matrix(angles)
% GET_ROTATION_MATRIX   determines a rotation matrix in euclidean space
%
% Syntax:
%   ROTATION_MATRIX = GET_ROTATION_MATRIX(ANGLES)
%
% Description:
%   ROTATION_MATRIX = GET_ROTATION_MATRIX(ANGLES) determines a 
%   ROTATION_MATRIX from three ANGLES (0 <= angle <= 360 degrees) with
%   respect to the x, y and z-axis respectively. ANGLES is a matrix with
%   the dimensions (1x3).
    if length(angles) == 3
        alpha = angles(1)*2*pi/360;
        beta  = angles(2)*2*pi/360;
        gamma = angles(3)*2*pi/360;
    
        R_x = [          1            0            0         0; ...
                         0   cos(alpha)  -sin(alpha)         0; ...
                         0   sin(alpha)   cos(alpha)         0; ...
                         0            0            0         1];
    
        R_y = [  cos(beta)            0    sin(beta)         0; ...
                         0            1            0         0; ...
                -sin(beta)            0    cos(beta)         0; ...
                         0            0            0         1];
    
        R_z = [ cos(gamma)  -sin(gamma)            0         0; ...
                sin(gamma)   cos(gamma)            0         0; ...
                         0            0            1         0; ...
                         0            0            0         1];
    
        rot_matrix = R_z * R_y * R_x;
    else    
        gamma = angles;
        
        R_z = [cosd(gamma) -sind(gamma)            0; ...
               sind(gamma)  cosd(gamma)            0; ...
                         0            0            1];
    
        rot_matrix = R_z;    
    end
end


function transformed_pointset = transform_pointset(pointset, transform)
    transformed_pointset = homo2cart((transform * cart2homo(pointset)')');
end


function transformed_surf = transform_surface(surface, transform)
    tmp                 = (transform * cart2homo(surface.vertices)')';
    transformed_surf    = struct('vertices', homo2cart(tmp), ...
                                 'faces', surface.faces);
end


%%
function [pc2, varargout] = rotatePointCloudAlongZ(pc, direction, varargin)
% This function rotates the point cloud of the highest eigenvector  
% to be along z-direction. 
% For the 2nd highest eigenvector, the user can choose to align 
% along x or y direction 

% parameters: 
% inputs: 
% pc: point cloud of size N * 3 size
% direction: 'x' or 'y' (align the 2nd highest eigenvector along 'x' or
% 'y' direction 
% output:
% pc2: rotated point cloud 
% 
% Based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z direction based on
% PCA, MATLAB Central File Exchange. Retrieved April 21, 2023.'

    if nargin > 2
        mean_pc = varargin{1};
        % Bring the point cloud center to the origin
        pc      = pc - mean_pc;
        u       = varargin{2};
    else
        mean_pc = mean(pc);
        % Bring the point cloud center to the origin
        pc      = pc - mean_pc;
        
        % Obtain the eigenvector of the highest eigenvalue
        u       = pcaEig(pc, 'max');
    end

    % This section align u normal vector along z-direction
       
    % Calculate the angles of the normal vector 
    [alpha, beta] = unitVectorToAngle(u);
    
    % Align the point cloud along x-axis followed by aligning along z-axis
    % YOU CAN REMOVE THE PI IF YOU WANT TO FLIP THE POINT CLOUD ALONG
    % Z-DIRECTION 
    [~, Ry, Rz] = rotationalMatrix(-alpha, pi-beta);
    pc2 = rotatePC(pc, Ry, Rz);
    
    % This section align v normal vector along x or y direction
    switch direction
        case 'x'
            offset = 0;
        case 'y'
            offset = pi/2;
    end
    
    if nargin > 2
        v = varargin{3};
    else
        % Obtain the eigenvector of the 2nd highest eigenvalue 
        v = pcaEig(pc2, 'middle');
    end

    % Calculate the angle of the projected v-vector along the xy-plane
    % with respect to the x-axis
    [alpha, ~] = unitVectorToAngle(v);
    
    % Calculate the rotational matrix for the angle
    [~, Ry, Rz] = rotationalMatrix(offset - alpha, 0);
    
    % Rotate the point cloud 
    pc2 = rotatePC(pc2, Ry, Rz);
    
    if nargout > 1
         varargout{1} = mean_pc;
         varargout{2} = u;
         varargout{3} = v;
    end
end


function [Rx, Ry, Rz] = rotationalMatrix(alpha, beta)
% Based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z direction based on
% PCA, MATLAB Central File Exchange. Retrieved April 21, 2023.'    
    Rx = [1 0 0; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)];
    Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
    Rz = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
end


function u = pcaEig(pc, magnitude)
% Based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z direction based on
% PCA, MATLAB Central File Exchange. Retrieved April 21, 2023.'
    
    % Obtain the covariance matrix
    covariance = cov([pc(:, 1) pc(:, 2) pc(:, 3)]);
   
    % Calculate the eigenvectors and obtain the normal vector
    [V, D] = eig(covariance);
    diagonalEigenvalues = diag(D);
    
    % Output the normal vector 
    
    % Sort the eigenvectors based on size of eigenvalues 
    [~, I] = sort(diagonalEigenvalues);
    V = V(:, I);
    
    switch magnitude
        case 'max'
            % Choose the eigenvector of the highest eigenvalue
            u = V(:, 3); 
        case 'middle'
            % Choose the eigenvector of the middle eigenvalue
            u = V(:, 2); 
        case 'min'
            % Choose the eigenvector of the lowest eigenvalue
            u = V(:, 1); 
    end    
end


function [alpha, beta] = unitVectorToAngle(u)
% Based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z direction based on
% PCA, MATLAB Central File Exchange. Retrieved April 21, 2023.'

    % Rotational angle between the projected u on the xy plane and the x-axis
    alpha = atan2(u(2), u(1)); 
    % Rotational angle between the u vector and the z-axis
    beta = atan2(sqrt(u(1)^2 + u(2)^2), u(3));
    
end


function pc2 = rotatePC(pc, Ry, Rz)
% Based on: 'Kin Sung Chan (2023). Align/Rotate Point Cloud Along Z direction based on
% PCA, MATLAB Central File Exchange. Retrieved April 21, 2023.'

    % Convert the point cloud to 3 * N format
    matrix = pc';
    % rotation around z axis to align point cloud along x axis 
    matrix2 = Rz*matrix;
    matrix2 = Ry*matrix2;
    % Ouput the point cloud 
    pc2 = matrix2';
end
