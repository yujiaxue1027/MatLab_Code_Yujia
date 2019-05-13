This is a preliminary version of Matlab code for column/row-wise "Decomposition of Space-Variant Blur in Image deconvolution".
Version: 0.8 (12.1.2016)
Author: Jan Kamenicky and Filip Sroubek, UTIA AV CR, v.v.i.
http://zoi.utia.cas.cz/decomposition

Description of m-files:
CWD algorithm is implemented in 'initCWD', 'iterU_CWD' and 'iterX_CWD'.
RWD algorithm is implemented in 'initRWD', 'iterU_RWD' and 'iterX_RWD'.
Usage is demonstrated in 'testWall' as an artificial experiment. At first, input images are degraded by SV convolution - cylinder blur with radius varying according to a linear gradient with different speed of change. Different SV deconvolution methods are then used to recover the original images for comparison (CWD, RWD, RWD using columns, RWD using bilinear interpolation instead of SVD).

The code will be improved and cleaned in the near future (February 2016).
