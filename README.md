# SRM-SIM
This is a Matlab implementation of spatially re-modulated structured illumination microscopy (SRM-SIM). It combines spatial re-modulation super-resoluion principle with Fourier domain filter and use measured illumination patterns, enabling high-speed and high-quality super-resolution SIM imaging. In addition, this strategy is also applicable to image reconstruction under distorted sinusoidal and other spatially uncorrelated illumination patterns. For detailed explanation, please see our paper [1]. To use this code, please also cite [1]. You can download the data used in the paper.

The method was proposed in the article "High-speed spatially re-modulated structured illumination microscopy" (https://doi.org/10.1364/xxxx)

The folder 'data' includes the raw simulated and experimental images.

The folder 'functions' includes the necessary subroutine;

'SRM_SIM.m' is the main algorithm to reconstruct a image in the CPU environment.

'SRM_SIM_GPU.m' is the main algorithm to reconstruct a image in the GPU environment.

'SRM_SIM_SpeckleMultifocal.m' is the main algorithm to reconstruct a image under speckle or multifocal illumination patters.

[1] Shijie Tu et al., "High-speed spatially remodulated structured illumination microscopy," xx.xx.xx
