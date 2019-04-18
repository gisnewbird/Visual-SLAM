% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 1415.630885504480100 ; 1418.157568236491900 ];

%-- Principal point:
cc = [ 647.590955842406740 ; 418.354233905402570 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.015451361846550 ; 0.111237201451761 ; 0.007216338486231 ; 0.002209335253680 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 27.694011958360033 ; 27.193799499804406 ];

%-- Principal point uncertainty:
cc_error = [ 19.431338856233072 ; 19.985715151281081 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.060633920125846 ; 0.503279336840889 ; 0.004523522772280 ; 0.005495260925043 ; 0.000000000000000 ];

%-- Image size:
nx = 1280;
ny = 720;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 10;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ 2.210372e+00 ; 2.103682e+00 ; -1.377674e-01 ];
Tc_1  = [ -7.628355e+01 ; -8.617164e+01 ; 4.131870e+02 ];
omc_error_1 = [ 1.081355e-02 ; 1.165621e-02 ; 2.324130e-02 ];
Tc_error_1  = [ 5.730500e+00 ; 5.802973e+00 ; 8.142799e+00 ];

%-- Image #2:
omc_2 = [ 1.997853e+00 ; 1.976899e+00 ; 1.269714e-01 ];
Tc_2  = [ -3.227494e+01 ; -8.452316e+01 ; 4.230671e+02 ];
omc_error_2 = [ 1.308400e-02 ; 1.132116e-02 ; 1.943718e-02 ];
Tc_error_2  = [ 5.847171e+00 ; 5.928776e+00 ; 8.613922e+00 ];

%-- Image #3:
omc_3 = [ -2.188458e+00 ; -2.036405e+00 ; 4.622760e-01 ];
Tc_3  = [ -1.281079e+02 ; -8.120017e+01 ; 4.603384e+02 ];
omc_error_3 = [ 1.380663e-02 ; 8.342058e-03 ; 2.372560e-02 ];
Tc_error_3  = [ 6.424435e+00 ; 6.566036e+00 ; 8.584268e+00 ];

%-- Image #4:
omc_4 = [ -2.153170e+00 ; -2.121084e+00 ; -1.945600e-03 ];
Tc_4  = [ -8.119246e+01 ; -8.956366e+01 ; 4.168205e+02 ];
omc_error_4 = [ 1.071340e-02 ; 1.050520e-02 ; 2.451230e-02 ];
Tc_error_4  = [ 5.795540e+00 ; 5.919917e+00 ; 8.142954e+00 ];

%-- Image #5:
omc_5 = [ 1.865381e+00 ; 1.827292e+00 ; -6.072784e-01 ];
Tc_5  = [ -7.725294e+01 ; -7.834414e+01 ; 4.746549e+02 ];
omc_error_5 = [ 9.719547e-03 ; 1.362851e-02 ; 1.838585e-02 ];
Tc_error_5  = [ 6.575141e+00 ; 6.696701e+00 ; 8.771149e+00 ];

%-- Image #6:
omc_6 = [ 1.763280e+00 ; 2.045492e+00 ; -1.723440e-01 ];
Tc_6  = [ 3.541968e+01 ; -9.383123e+01 ; 6.700832e+02 ];
omc_error_6 = [ 1.186866e-02 ; 1.332924e-02 ; 2.008539e-02 ];
Tc_error_6  = [ 9.261169e+00 ; 9.417577e+00 ; 1.340228e+01 ];

%-- Image #7:
omc_7 = [ 2.094634e+00 ; 2.012057e+00 ; 1.949794e-01 ];
Tc_7  = [ -3.355551e+00 ; -3.539211e+01 ; 6.063796e+02 ];
omc_error_7 = [ 1.500982e-02 ; 1.109124e-02 ; 2.295356e-02 ];
Tc_error_7  = [ 8.317950e+00 ; 8.522138e+00 ; 1.229017e+01 ];

%-- Image #8:
omc_8 = [ -1.882702e+00 ; -2.210035e+00 ; 4.383333e-01 ];
Tc_8  = [ -9.101592e-01 ; -7.997921e+01 ; 6.264811e+02 ];
omc_error_8 = [ 1.189182e-02 ; 1.234174e-02 ; 2.576158e-02 ];
Tc_error_8  = [ 8.623966e+00 ; 8.794925e+00 ; 1.188207e+01 ];

%-- Image #9:
omc_9 = [ 2.109612e+00 ; 1.736796e+00 ; -7.977224e-01 ];
Tc_9  = [ -9.323661e+01 ; -7.077717e+01 ; 6.461867e+02 ];
omc_error_9 = [ 9.923244e-03 ; 1.403531e-02 ; 2.052114e-02 ];
Tc_error_9  = [ 8.945650e+00 ; 9.137428e+00 ; 1.184833e+01 ];

%-- Image #10:
omc_10 = [ 1.853428e+00 ; 1.708387e+00 ; -6.873696e-01 ];
Tc_10  = [ -6.048933e+01 ; -1.219389e+02 ; 6.361252e+02 ];
omc_error_10 = [ 1.077956e-02 ; 1.447083e-02 ; 1.826023e-02 ];
Tc_error_10  = [ 8.867928e+00 ; 8.989887e+00 ; 1.197554e+01 ];

