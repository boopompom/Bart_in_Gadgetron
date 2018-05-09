
# Integrate Bart into Gadgetron
1. BartGccGadget implements Geometric Coil Compression (GCC) [1] for Cartesian 3D data. 
2. BartReconGadget calls ESPIRiT calibration and PICS reconstruciton provided in Bart to implement L1-ESPIRiT reconstruction of Cartesian 3D data. The communication between Bart and Gadgetron is through a user-defined shell script file (which can be used/tested without Gadgetron). The data write/read is implemented via .cfl/.hdr files.


Problem: ecalib commmand run slowly in Gadgetron for large 3D datasets.

Reference
[1] Tao Zhang, JM Pauly, SS Vasanawala, Michael Lustig. Coil Compression for Accelerated Imaging with Cartesian Sampling. Magn Reson Med, 2013, 69:571-582.
[2] Martin Uecker, Peng Lai, MJ Murphy, Patrick Virtue, Micahel Elad, JM Pauly, SS Vasanawala, Michael Lustig. ESPIRiT-An Eigenvalue Approach to Autocalibrating Parallel MRI: where SENSE meets GRAPPA. Magn Reson Med, 2014, 71:990-1001.
[3] http://gadgetron.github.io/
[4] http://mrirecon.github.io/bart/
