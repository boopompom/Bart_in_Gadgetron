# Integrate Bart into Gadgetron
1. BartGccGadget implements Geometric Coil Compression (GCC) [1] for Cartesian 3D data. 
2. BartReconGadget for L1-ESPIRiT Reconstruction of Cartesian 3D data. The communication between Bart and Gadgetron is through a user-defined script file. The data write/read is implemented via .cfl/.hdr files.

Reference
[1] Tao Zhang, JM Pauly, SS Vasanawala, Michael Lustig. Coil Compression for Accelerated Imaging with Cartesian Sampling. Magn Reson Med, 2013, 69:571-582.
[2] Martin Uecker, Peng Lai, MJ Murphy, Patrick Virtue, Micahel Elad, JM Pauly, SS Vasanawala, Michael Lustig. ESPIRiT-An Eigenvalue Approach to Autocalibrating Parallel MRI: where SENSE meets GRAPPA. Magn Reson Med, 2014, 71:990-1001.
[3] http://gadgetron.github.io/
[4] http://mrirecon.github.io/bart/
