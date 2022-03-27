% Rafael Villamor Lora 
% May 5, 2020 (Last updated 7/01/2020) 
% FROM 2D MAPS TO STL 
close all
%%-----------------------------------------------------------------------%% 

%            GENERATE A JOINT MODEL USING BROWN'S (1995) MODEL            % 

%%-----------------------------------------------------------------------%% 

% INPUTS 
L         = 25.4*3;     % Fracture length 
H         = 0.8;        % Hurst exponent. Determines fractal dimension D=3-H                  [0.45<H<0.85] 
roughness = 0.005*L;    % root-mean-square roughness (~ standard deviation of heights)        [same units as L] 
mismatch  = 0.05;       % Mismatch length scale (wavelength) as a fraction of fracture size   [0<Mismatch<1] 
N         = 2^10;       % Number of divisions along fracture side                             [power of 2] 
aniso     = 1.0;        % anisotropy ratio                                                    [0<Aniso<1] 
seed      = 2;          % Random number seed (phase generation)                               [+ve integer] 
dx        = L/N;        % Grid spacing                                                        [same units as L] 
lambda_0  = 0.25;       % Roll-ff wavelenth                                                   [0<lambda_0<1] 
model     = 'bilinear'; % Power spectal density model                                         ['linear', 'linear', 'bilinear', 'smooth'] 

  

% GENERATE FRACTURE 

[Z{2}, Z{1}, Zap] = RSG_brown1995(H, roughness, mismatch, N, aniso, seed, lambda_0, model); 

  

% TRANSLATE SURFACES UNTIL THEY ARE IN CONTACT 

Z{1}      = Z{1} - min(Z{1}(:)); 
shift_Z2  = Z{2} - Z{1}; 
Z{2}      = Z{2} - min(shift_Z2(:)); 

  

% TRANSLATE SURFACE SO THE APERTURE VOLUME IT'S LOCATED AT Z = 0 

shift_Z12 = (mean(Z{1}(:)) + mean(Z{2}(:)) ) / 2; 

Z{1}      = Z{1} - shift_Z12; 

Z{2}      = Z{2} - shift_Z12; 

  

% APERTURE FIELD 

Z{3}      = Zap; 

  

%%-----------------------------------------------------------------------%% 

%                             GENERATE BODIES                             % 

%%-----------------------------------------------------------------------%% 

  

% This section uses the surf2solid() and STLwrite() functions from Sven: 

% https://www.mathworks.com/matlabcentral/fileexchange/42876-surf2solid-make-a-solid-volume-from-a-surface-for-3d-printing 

% https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-write-ascii-or-binary-stl-files 

% 

% See some examples of MATLAB 3D printing @ 

% Reference: https://efcms.engr.utk.edu/ef230-2020-06/modules/3dprinting/matlab3d.php 

  

% STL PARAMTERES 

geometry = 'prism';                                                        % Specimen geometry: 'prism' or 'cylinder' 

Nstl     = N;                                                              % Resolution of the STL file [integer <= N] 

stlFiles = {'Results/lowerFracture.stl',  'Results/upperFracture.stl';...  % File names 

            'Results/apertureVolume.stl', 'Results/apertureMold.stl'};     

  

% SPECIMEN DIMENSIONS 

switch geometry 

    case 'prism' 

        Spcm_width    = dx * N;                                            % Specimen Width  [same units as L] 

        Spcm_depth    = dx * N/2;                                          % Specimen Depth  [same units as L] 

        Spcm_height   = 24.5*0.25;                                         % Specimen height [same units as L] 

    case 'cylinder' 

        Spcm_diameter = dx * N;                                            % Specimen diameter [same units as L] 

        Spcm_height   = dx * 2 * N;                                        % Specimen height   [same units as L] 

    otherwise 

        error('Geometry not supported'); 

end 

  

for i = 1:2 % Repeat loop for bottom (i = 1) and top (i = 2) 

    % CHANGE THE RESOLUTION OF THE STL (in case the file is too heavy) 

    fracture = Z{i};                                                       % Fracture geometry 

    fracture = imresize(fracture, [Nstl Nstl]);                            % Resize 

     

    % CROP TO DESIRED DIMENSIONS 

    switch geometry 

        case 'prism' 

            fracture_stl = fracture(1:Spcm_width/dx, 1:Spcm_depth/dx);                                                              % Region of interest 

            [XX, YY]     = meshgrid(linspace(0, Spcm_depth, size(fracture_stl,2)), linspace(0, Spcm_width, size(fracture_stl,1)));  % X-Y matrices 

        case 'cylinder' 

            fracture_stl = fracture;                                                                                                % Region of interest 

            [XX, YY]     = meshgrid(linspace(0, Spcm_diameter, Nstl), linspace(0, Spcm_diameter, Nstl));                            % X-Y matrices   

        otherwise 

            error('Geometry not supported'); 

    end 

     

    % SHIFT THE FRACTURE TO THE APROPRIATE ELEVATION 

    midPoint    = Spcm_height / 2; 

    if     i == 1 % Bottom surface 

        ZZ = fracture_stl + midPoint;                                      % Translate 

    elseif i == 2 % Top surface 

        ZZ = fracture_stl;                                                  

        R  = roty(180) * [XX(:) YY(:) ZZ(:)]';                             % Flip the top surface 

        XX = reshape(R(1,:), size(XX,1), size(XX,2));                      % New x-coordinate matrix 

        YY = reshape(R(2,:), size(YY,1), size(YY,2));                      % New y-coordinate matrix 

        ZZ = reshape(R(3,:), size(ZZ,1), size(ZZ,2)) + midPoint;           % New z-coordinate matrix + Translate 

    end 

     

    % CROP TO DESIRED GEOMETRY 

    switch geometry 

        case 'prism' 

        case 'cylinder' 

            [~,rho]     = cart2pol(XX - mean(XX(1,:)),YY - mean(YY(:,1))); % Convert to polar coordinates 

            outside     = find(rho > Spcm_diameter/2);                     % Points outside the cylinder 

            ZZ(outside) = 0;                                               % Set points outside = 0 

        otherwise 

            error('Geometry not supported'); 

    end 

         

    % CONVERT SURFACE TO SOLID AND SAVE STL FILE 

FV{i} = surf2solid(XX, YY, ZZ, 'elevation', 0);                        %#ok<SAGROW> % Convert to a closed solid 
stlwrite(stlFiles{i}, FV{i});                                          % Save as .STL 
    % SAVE SURFACES FOR LATER 
    if     i == 1 % Bottom surface 

        surfaces{1} = ZZ - midPoint;                                       % Move back to the origin 

    elseif i == 2 % Top surface 

        ZZ = ZZ - midPoint;                                                % Move back to the origin 

        R  = roty(180) * [XX(:) YY(:) ZZ(:)]';                             % Flip the top surface 

        XX = reshape(R(1,:), size(XX,1), size(XX,2));                

        YY = reshape(R(2,:), size(YY,1), size(YY,2)); 

        ZZ = reshape(R(3,:), size(ZZ,1), size(ZZ,2)); 

        surfaces{2} = ZZ; 

    end 

     

end 

  

% APERTURE VOLUME (in case one wants to do 3D simulations) 

FV{3} = surf2solid(XX, YY, surfaces{2}, 'elevation', surfaces{1});         % Generate the solid volume 

stlwrite(stlFiles{3}, FV{3});                                              % Save as .STL 

  

% APERTURE MOLD 

apertureMold = - (surfaces{2} - surfaces{1});                              % The mold is the negative of the aperture field 

apertureMold = apertureMold - min(apertureMold(:));                        % Translate so the lowest point is at the base 

FV{4}        = surf2solid(XX, YY, apertureMold, 'elevation', 0);           % Generate the solid volume 

stlwrite(stlFiles{4}, FV{4});                                              % Save as .STL 

  

%%-----------------------------------------------------------------------%% 

%                               PLOT RESULTS                              % 

%%-----------------------------------------------------------------------%% 

  

% PLOT JOINT 

figure('Position', [1922 41 577 953]) 

% Lower surface 

trisurf(FV{1}.faces, FV{1}.vertices(:,1), FV{1}.vertices(:,2),  FV{1}.vertices(:,3)); hold on   % Plot lower volume 

% Upper surface 

R = roty(180) * [FV{2}.vertices(:,1) FV{2}.vertices(:,2) FV{2}.vertices(:,3)]';                 % Flip the upper body  

trisurf(FV{2}.faces, R(1,:)', R(2,:)',  R(3,:)' + 2*midPoint);                                  % Plot upper volume 

axis equal vis3d, shading interp, camlight left, lighting gouraud, colormap white               % Make everything look pretty  

switch geometry                                                                                 % Set the scale 

    case 'prism' 

        axis([0 Spcm_depth    0 Spcm_width    0 Spcm_height]); 

    case 'cylinder' 

        axis([0 Spcm_diameter 0 Spcm_diameter 0 Spcm_height]); 

end 

  

% PLOT APERTURE VOLUME 

figure('Position', [2501 41 577 953]) 

trisurf(FV{3}.faces, FV{3}.vertices(:,1), FV{3}.vertices(:,2), FV{3}.vertices(:,3) + midPoint); % Plot aperture volume 

axis equal vis3d, shading interp, camlight left, lighting gouraud, colormap white               % Make everything look pretty  

switch geometry                                                                                 % Set the scale 

    case 'prism' 

        axis([0 Spcm_depth    0 Spcm_width    0 Spcm_height]); 

    case 'cylinder' 

        axis([0 Spcm_diameter 0 Spcm_diameter 0 Spcm_height]); 

end 

  

% PLOT APERTURE MOLD 

figure('Position', [3079 41 577 953]) 

trisurf(FV{4}.faces, FV{4}.vertices(:,1), FV{4}.vertices(:,2), FV{4}.vertices(:,3)); % Plot aperture volume 

axis equal vis3d, shading interp, camlight left, lighting gouraud, colormap white    % Make everything look pretty  

  

%%-----------------------------------------------------------------------%% 

%                               SAVE RESULTS                              % 

%%-----------------------------------------------------------------------%% 

generatedSurface.L         = L;          % Fracture length 

generatedSurface.H         = H;          % Hurst exponent. Determines fractal dimension D=3-H                   

generatedSurface.roughness = roughness;  % root-mean-square roughness (i.e.variance of heights)                 

generatedSurface.mismatch  = mismatch;   % Mismatch length scale (wavelength) as a fraction of fracture size    

generatedSurface.N         = N;          % Number of divisions along fracture side                             

generatedSurface.aniso     = aniso;      % anisotropy ratio                                                     

generatedSurface.seed      = seed;       % Random number seed (phase generation)                                

generatedSurface.dx        = dx;         % Grid spacing                                                         

generatedSurface.lambda_0  = lambda_0;   % Roll-ff wavelenth                                                    

generatedSurface.model     = model;      % Power spectal density model                                          

generatedSurface.seed      = seed;       % Random number seed (phase generation)                                

generatedSurface.Z         = Z;          % Generated surfaces 

  

generatedSurface.geometry  = geometry;   % Specimen geometry: 'prism' or 'cylinder' 

generatedSurface.Nstl      = Nstl;       % Resolution of the STL file [integer <= N] 

switch geometry 

    case 'prism' 

        generatedSurface.Spcm_width    = Spcm_width;     % Specimen Width   

        generatedSurface.Spcm_depth    = Spcm_depth;     % Specimen Depth   

        generatedSurface.Spcm_height   = Spcm_height;    % Specimen height  

    case 'cylinder' 

        generatedSurface.Spcm_diameter = Spcm_diameter;  % Specimen diameter  

        generatedSurface.Spcm_height   = Spcm_height;    % Specimen height  

    otherwise 

        error('Geometry not supported'); 

end 

generatedSurface.surfaces  = surfaces;   % Lower surface, upper surface 

  

save('Results/generatedSurface.mat', 'generatedSurface'); 

 