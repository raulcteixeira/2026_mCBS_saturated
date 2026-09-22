Attenuation_Raul.txt: the OD created by the cloud, taking into account saturation, for input and reflected 
beam, for each s_0.

atomic_movement_...: no absorption. Corresponds to equal power in incoming and outcoming, that is, 
    the atoms at the end of the cloud.
atomic_movement_abs_1_... : corresponds to points in the beginning of the cloud, that is, with a 
    difference of power for both beams that correspond to absorption in incoming + reflected beam
atomic_movement_abs_2_... : corresponds to points in the middle of the cloud, that is, with a 
    difference of power for both beams that correspond to half the absorption in incoming + reflected beam

The three file groups above are used for the paper. The timevector is from 0 to 10us in steps of 0.5 us

Below, all the same, but with timevector is from 0 to 2us in steps of 0.1 us
small_time_scale/atomic_movement_..., atomic_movement_abs_1..., atomic_movement_abs_2... 