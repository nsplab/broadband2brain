This code compares the performance of
1) Fast Nonnegative Deconvolution (Vogelstein et al.)
2) Iterative ML
applied to the problem of reconstructing spike times from calcium imaging

Overview of algorithms:

1) Outputs vector n_t where each entries represents probability of spike
in that time bin. All spikes are assumed to have the same size. Need to
apply threshold afterwards if you want a spike train.

2) Assumes there are exactly K spikes, finds approximate maximum
likelihood estimate for their times and amplitudes (t, c). Small spikes
can be rejected using threshold.








Imaging-SNR-Data:

  Variables:
    ExpDat20080716a24    329215    1960
    ExpDat20080716a26    311607    1770
    ExpDat20080716a28    272817    1653
    ExpDat20080716a31    541878    3415
    ExpDat20080716a35    378698    2316
    ExpDat20080716a39    371883     829
    ExpDat20080721a10    444245    2256
    ExpDat20080721a12    543249    2930
    ExpDat20080721a7     676197    2674       *
    ExpDat20080721a8     456437    2411       *
    ExpDat20080801b4     483399       1
    ExpDat20080801b5     641599       1
                            A        B    Has FOOPSI


  Each of the variables is a structure with the fields:
    chanDev1_ai0_VoltageCh1: [Ax1 single]
    chanDev1_ai1_CurrentCh1: [Ax1 single]
    chanDev1_ai2_CameraSync: [Ax1 single]
    time: [Ax1 double]
    Fluorescence: [1xB double]
    FOOPSI: [1xB double]
    FluorescenceTime: [Bx1 double]
    Analysis: [1x1 struct]
