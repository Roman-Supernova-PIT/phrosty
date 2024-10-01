import cupy as cp
import numpy as np
from sfft.utils.SFFTSolutionReader import Realize_MatchingKernel
from sfft.PureCupyCustomizedPacket import PureCupy_Customized_Packet
from SpaceSFFTPurePack.sfftutils.PatternRotationCalculator import PatternRotation_Calculator
from SpaceSFFTPurePack.sfftutils.PureCupyFFTKits import PureCupy_FFTKits
from SpaceSFFTPurePack.sfftutils.PureCupyDeCorrelationCalculator import PureCupy_DeCorrelation_Calculator

# loacal imports
from SpaceSFFTPurePack.ResampKits import Cupy_ZoomRotate
from SpaceSFFTPurePack.ResampKits import Cupy_Resampling

__last_update__ = "2024-09-22"
__author__ = "Lei Hu <leihu@andrew.cmu.edu>"

class SpaceSFFT_CupyFlow:
    @staticmethod
    def SSCF(hdr_REF, hdr_oSCI, PixA_REF_GPU, PixA_oSCI_GPU, PixA_REF_DMASK_GPU, PixA_oSCI_DMASK_GPU, 
        PSF_REF_GPU, PSF_oSCI_GPU, GKerHW=9, KerPolyOrder=2, BGPolyOrder=0, ConstPhotRatio=True, 
        CUDA_DEVICE_4SUBTRACT='0', GAIN=1.0):
        """Run A Cupy WorkFlow for SFFT subtraction"""

        assert PixA_REF_GPU.flags['C_CONTIGUOUS']
        assert PixA_oSCI_GPU.flags['C_CONTIGUOUS']
        
        assert PixA_REF_DMASK_GPU.flags['C_CONTIGUOUS']
        assert PixA_oSCI_DMASK_GPU.flags['C_CONTIGUOUS']
        
        assert PSF_REF_GPU.flags['C_CONTIGUOUS']
        assert PSF_oSCI_GPU.flags['C_CONTIGUOUS']

        assert "BKG_SIG" in hdr_REF
        assert "BKG_SIG" in hdr_oSCI

        # * step 0. run resampling for input science image, mask and PSF
        CR = Cupy_Resampling(RESAMP_METHOD="BILINEAR", VERBOSE_LEVEL=1)
        XX_proj_GPU, YY_proj_GPU = CR.resamp_projection_sip(hdr_obj=hdr_oSCI, hdr_targ=hdr_REF, NSAMP=1024, RANDOM_SEED=10086)

        PixA_Eobj_GPU, EProjDict = CR.frame_extension(XX_proj_GPU=XX_proj_GPU, YY_proj_GPU=YY_proj_GPU, 
            PixA_obj_GPU=PixA_oSCI_GPU, PAD_FILL_VALUE=0., NAN_FILL_VALUE=0.)
        PixA_SCI_GPU = CR.resampling(PixA_Eobj_GPU=PixA_Eobj_GPU, EProjDict=EProjDict, USE_SHARED_MEMORY=False)

        PixA_Eobj_GPU, EProjDict = CR.frame_extension(XX_proj_GPU=XX_proj_GPU, YY_proj_GPU=YY_proj_GPU, 
            PixA_obj_GPU=PixA_oSCI_DMASK_GPU, PAD_FILL_VALUE=0., NAN_FILL_VALUE=0.)
        PixA_SCI_DMASK_GPU = CR.resampling(PixA_Eobj_GPU=PixA_Eobj_GPU, EProjDict=EProjDict, USE_SHARED_MEMORY=False)
        BlankMask_GPU = PixA_SCI_GPU == 0.

        PATTERN_ROTATE_ANGLE = PatternRotation_Calculator.PRC(hdr_obj=hdr_oSCI, hdr_targ=hdr_REF)
        PSF_SCI_GPU = Cupy_ZoomRotate.CZR(PixA_obj_GPU=PSF_oSCI_GPU, ZOOM_SCALE_X=1., ZOOM_SCALE_Y=1., \
            OUTSIZE_PARIRY_X='UNCHANGED', OUTSIZE_PARIRY_Y='UNCHANGED', PATTERN_ROTATE_ANGLE=PATTERN_ROTATE_ANGLE, \
            RESAMP_METHOD='BILINEAR', PAD_FILL_VALUE=0., NAN_FILL_VALUE=0., THREAD_PER_BLOCK=8, \
            USE_SHARED_MEMORY=False, VERBOSE_LEVEL=2)

        # * step 1. cross convolution
        PixA_CREF_GPU = PureCupy_FFTKits.FFT_CONVOLVE(PixA_Inp_GPU=PixA_REF_GPU, KERNEL_GPU=PSF_SCI_GPU, 
            PAD_FILL_VALUE=0., NAN_FILL_VALUE=None, NORMALIZE_KERNEL=True, FORCE_OUTPUT_C_CONTIGUOUS=True, FFT_BACKEND="Cupy")

        PixA_CSCI_GPU = PureCupy_FFTKits.FFT_CONVOLVE(PixA_Inp_GPU=PixA_SCI_GPU, KERNEL_GPU=PSF_REF_GPU,
            PAD_FILL_VALUE=0., NAN_FILL_VALUE=None, NORMALIZE_KERNEL=True, FORCE_OUTPUT_C_CONTIGUOUS=True, FFT_BACKEND="Cupy")

        # * step 2. sfft subtraction
        LYMASK_BKG_GPU = cp.logical_or(PixA_REF_DMASK_GPU == 0, PixA_SCI_DMASK_GPU < 0.1)   # background-mask

        NaNmask_CREF_GPU = cp.isnan(PixA_CREF_GPU)
        NaNmask_CSCI_GPU = cp.isnan(PixA_CSCI_GPU)
        if NaNmask_CREF_GPU.any() or NaNmask_CSCI_GPU.any():
            NaNmask_GPU = cp.logical_or(NaNmask_CREF_GPU, NaNmask_CSCI_GPU)
            ZeroMask_GPU = cp.logical_or(NaNmask_GPU, LYMASK_BKG_GPU)
        else:
            ZeroMask_GPU = LYMASK_BKG_GPU

        PixA_mCREF_GPU = PixA_CREF_GPU.copy()
        PixA_mCREF_GPU[ZeroMask_GPU] = 0.

        PixA_mCSCI_GPU = PixA_CSCI_GPU.copy()
        PixA_mCSCI_GPU[ZeroMask_GPU] = 0.

        # trigger sfft subtraction
        ForceConv = 'REF'    # here force to REF
        Solution_GPU, PixA_DIFF_GPU = PureCupy_Customized_Packet.PCCP(
            PixA_REF_GPU=PixA_CREF_GPU, PixA_SCI_GPU=PixA_CSCI_GPU, PixA_mREF_GPU=PixA_mCREF_GPU, PixA_mSCI_GPU=PixA_mCSCI_GPU, 
            ForceConv=ForceConv, GKerHW=GKerHW, KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, 
            CUDA_DEVICE_4SUBTRACT=CUDA_DEVICE_4SUBTRACT
        )
        PixA_DIFF_GPU[BlankMask_GPU] = 0.

        # * step 3. perform decorrelation in Fourier domain
        # extract matching kernel at the center
        N0, N1 = PixA_DIFF_GPU.shape
        L0, L1 = 2*GKerHW + 1, 2*GKerHW + 1
        DK = KerPolyOrder
        Fpq = int((BGPolyOrder+1)*(BGPolyOrder+2)/2)
        XY_q = np.array([[N0/2.+0.5, N1/2.+0.5]])

        Solution = cp.asnumpy(Solution_GPU)
        MATCH_KERNEL_GPU = cp.array(Realize_MatchingKernel(XY_q=XY_q).FromArray(
            Solution=Solution, N0=N0, N1=N1, L0=L0, L1=L1, DK=DK, Fpq=Fpq
        )[0], dtype=cp.float64)

        # do decorrelation
        # Note: we ignored the background noise change by resampling
        BKGSIG_SCI = hdr_oSCI['BKG_SIG']  
        BKGSIG_REF = hdr_REF['BKG_SIG']

        FKDECO_GPU = PureCupy_DeCorrelation_Calculator.PCDC(NX_IMG=N0, NY_IMG=N1, KERNEL_GPU_JQueue=[PSF_REF_GPU], 
            BKGSIG_JQueue=[BKGSIG_SCI], KERNEL_GPU_IQueue=[PSF_SCI_GPU], BKGSIG_IQueue=[BKGSIG_REF], 
            MATCH_KERNEL_GPU=MATCH_KERNEL_GPU, REAL_OUTPUT=False, REAL_OUTPUT_SIZE=None, 
            NORMALIZE_OUTPUT=True, VERBOSE_LEVEL=2)

        FPixA_DIFF_GPU = cp.fft.fft2(PixA_DIFF_GPU)
        PixA_DCDIFF_GPU = cp.fft.ifft2(FPixA_DIFF_GPU * FKDECO_GPU).real

        # * step 4. noise decorrelation & SNR estimation 
        # roughly estimate the SNR map for the decorrelated difference image
        # WARNING: the noise propagation is highly simplified.

        PixA_varREF_GPU = cp.clip(PixA_REF_GPU/GAIN, a_min=0.0, a_max=None) + BKGSIG_REF**2
        PixA_varSCI_GPU = cp.clip(PixA_SCI_GPU/GAIN, a_min=0.0, a_max=None) + BKGSIG_SCI**2
        PixA_NDIFF_GPU = cp.sqrt(PixA_varREF_GPU + PixA_varSCI_GPU)
        PixA_DSNR_GPU = PixA_DCDIFF_GPU / PixA_NDIFF_GPU

        return PixA_DIFF_GPU, PixA_DCDIFF_GPU, PixA_DSNR_GPU

class SpaceSFFT_CupyFlow_NVTX:
    @staticmethod
    def SSCFN(hdr_REF, hdr_oSCI, PixA_REF_GPU, PixA_oSCI_GPU, PixA_REF_DMASK_GPU, PixA_oSCI_DMASK_GPU, 
        PSF_REF_GPU, PSF_oSCI_GPU, GKerHW=9, KerPolyOrder=2, BGPolyOrder=0, ConstPhotRatio=True, 
        CUDA_DEVICE_4SUBTRACT='0', GAIN=1.0):
        """Run A Cupy WorkFlow for SFFT subtraction"""

        assert PixA_REF_GPU.flags['C_CONTIGUOUS']
        assert PixA_oSCI_GPU.flags['C_CONTIGUOUS']
        
        assert PixA_REF_DMASK_GPU.flags['C_CONTIGUOUS']
        assert PixA_oSCI_DMASK_GPU.flags['C_CONTIGUOUS']
        
        assert PSF_REF_GPU.flags['C_CONTIGUOUS']
        assert PSF_oSCI_GPU.flags['C_CONTIGUOUS']

        assert "BKG_SIG" in hdr_REF
        assert "BKG_SIG" in hdr_oSCI

        import nvtx
        
        # * step 0. run resampling for input science image, mask and PSF
        with nvtx.annotate("0-resamp", color="#FFD900"):
            CR = Cupy_Resampling(RESAMP_METHOD="BILINEAR", VERBOSE_LEVEL=1)
            XX_proj_GPU, YY_proj_GPU = CR.resamp_projection_sip(hdr_obj=hdr_oSCI, hdr_targ=hdr_REF, NSAMP=1024, RANDOM_SEED=10086)

            PixA_Eobj_GPU, EProjDict = CR.frame_extension(XX_proj_GPU=XX_proj_GPU, YY_proj_GPU=YY_proj_GPU, 
                PixA_obj_GPU=PixA_oSCI_GPU, PAD_FILL_VALUE=0., NAN_FILL_VALUE=0.)
            PixA_SCI_GPU = CR.resampling(PixA_Eobj_GPU=PixA_Eobj_GPU, EProjDict=EProjDict, USE_SHARED_MEMORY=False)

            PixA_Eobj_GPU, EProjDict = CR.frame_extension(XX_proj_GPU=XX_proj_GPU, YY_proj_GPU=YY_proj_GPU, 
                PixA_obj_GPU=PixA_oSCI_DMASK_GPU, PAD_FILL_VALUE=0., NAN_FILL_VALUE=0.)
            PixA_SCI_DMASK_GPU = CR.resampling(PixA_Eobj_GPU=PixA_Eobj_GPU, EProjDict=EProjDict, USE_SHARED_MEMORY=False)
            BlankMask_GPU = PixA_SCI_GPU == 0.

            PATTERN_ROTATE_ANGLE = PatternRotation_Calculator.PRC(hdr_obj=hdr_oSCI, hdr_targ=hdr_REF)
            PSF_SCI_GPU = Cupy_ZoomRotate.CZR(PixA_obj_GPU=PSF_oSCI_GPU, ZOOM_SCALE_X=1., ZOOM_SCALE_Y=1., \
                OUTSIZE_PARIRY_X='UNCHANGED', OUTSIZE_PARIRY_Y='UNCHANGED', PATTERN_ROTATE_ANGLE=PATTERN_ROTATE_ANGLE, \
                RESAMP_METHOD='BILINEAR', PAD_FILL_VALUE=0., NAN_FILL_VALUE=0., THREAD_PER_BLOCK=8, \
                USE_SHARED_MEMORY=False, VERBOSE_LEVEL=2)
        
        # * step 1. cross convolution
        with nvtx.annotate("1-crossConv", color="#C51B8A"):
            PixA_CREF_GPU = PureCupy_FFTKits.FFT_CONVOLVE(PixA_Inp_GPU=PixA_REF_GPU, KERNEL_GPU=PSF_SCI_GPU, 
                PAD_FILL_VALUE=0., NAN_FILL_VALUE=None, NORMALIZE_KERNEL=True, FORCE_OUTPUT_C_CONTIGUOUS=True, FFT_BACKEND="Cupy")

            PixA_CSCI_GPU = PureCupy_FFTKits.FFT_CONVOLVE(PixA_Inp_GPU=PixA_SCI_GPU, KERNEL_GPU=PSF_REF_GPU,
                PAD_FILL_VALUE=0., NAN_FILL_VALUE=None, NORMALIZE_KERNEL=True, FORCE_OUTPUT_C_CONTIGUOUS=True, FFT_BACKEND="Cupy")

        # * step 2. sfft subtraction
        with nvtx.annotate("2-sfft", color="#233ED9"):
            LYMASK_BKG_GPU = cp.logical_or(PixA_REF_DMASK_GPU == 0, PixA_SCI_DMASK_GPU < 0.1)   # background-mask

            NaNmask_CREF_GPU = cp.isnan(PixA_CREF_GPU)
            NaNmask_CSCI_GPU = cp.isnan(PixA_CSCI_GPU)
            if NaNmask_CREF_GPU.any() or NaNmask_CSCI_GPU.any():
                NaNmask_GPU = cp.logical_or(NaNmask_CREF_GPU, NaNmask_CSCI_GPU)
                ZeroMask_GPU = cp.logical_or(NaNmask_GPU, LYMASK_BKG_GPU)
            else:
                ZeroMask_GPU = LYMASK_BKG_GPU

            PixA_mCREF_GPU = PixA_CREF_GPU.copy()
            PixA_mCREF_GPU[ZeroMask_GPU] = 0.

            PixA_mCSCI_GPU = PixA_CSCI_GPU.copy()
            PixA_mCSCI_GPU[ZeroMask_GPU] = 0.

            # trigger sfft subtraction
            ForceConv = 'REF'  # force to REF
            Solution_GPU, PixA_DIFF_GPU = PureCupy_Customized_Packet.PCCP(
                PixA_REF_GPU=PixA_CREF_GPU, PixA_SCI_GPU=PixA_CSCI_GPU, PixA_mREF_GPU=PixA_mCREF_GPU, PixA_mSCI_GPU=PixA_mCSCI_GPU, 
                ForceConv=ForceConv, GKerHW=GKerHW, KerPolyOrder=KerPolyOrder, BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, 
                CUDA_DEVICE_4SUBTRACT=CUDA_DEVICE_4SUBTRACT
            )
            PixA_DIFF_GPU[BlankMask_GPU] = 0.

        # * step 3. perform decorrelation in Fourier domain
        with nvtx.annotate("3-decorr", color="#65BD63"):
            # extract matching kernel at the center
            N0, N1 = PixA_DIFF_GPU.shape
            L0, L1 = 2*GKerHW + 1, 2*GKerHW + 1
            DK = KerPolyOrder
            Fpq = int((BGPolyOrder+1)*(BGPolyOrder+2)/2)
            XY_q = np.array([[N0/2.+0.5, N1/2.+0.5]])

            Solution = cp.asnumpy(Solution_GPU)
            MATCH_KERNEL_GPU = cp.array(Realize_MatchingKernel(XY_q=XY_q).FromArray(
                Solution=Solution, N0=N0, N1=N1, L0=L0, L1=L1, DK=DK, Fpq=Fpq
            )[0], dtype=cp.float64)

            # do decorrelation
            # Note: we ignored the background noise change by resampling
            BKGSIG_SCI = hdr_oSCI['BKG_SIG']  
            BKGSIG_REF = hdr_REF['BKG_SIG']

            FKDECO_GPU = PureCupy_DeCorrelation_Calculator.PCDC(NX_IMG=N0, NY_IMG=N1, KERNEL_GPU_JQueue=[PSF_REF_GPU], 
                BKGSIG_JQueue=[BKGSIG_SCI], KERNEL_GPU_IQueue=[PSF_SCI_GPU], BKGSIG_IQueue=[BKGSIG_REF], 
                MATCH_KERNEL_GPU=MATCH_KERNEL_GPU, REAL_OUTPUT=False, REAL_OUTPUT_SIZE=None, 
                NORMALIZE_OUTPUT=True, VERBOSE_LEVEL=2)

            FPixA_DIFF_GPU = cp.fft.fft2(PixA_DIFF_GPU)
            PixA_DCDIFF_GPU = cp.fft.ifft2(FPixA_DIFF_GPU * FKDECO_GPU).real
        
        # * step 4. noise decorrelation & SNR estimation 
        with nvtx.annotate("4-dsnr", color="#00D0FF"):
            # roughly estimate the SNR map for the decorrelated difference image
            # WARNING: the noise propagation is highly simplified.

            PixA_varREF_GPU = cp.clip(PixA_REF_GPU/GAIN, a_min=0.0, a_max=None) + BKGSIG_REF**2
            PixA_varSCI_GPU = cp.clip(PixA_SCI_GPU/GAIN, a_min=0.0, a_max=None) + BKGSIG_SCI**2
            PixA_NDIFF_GPU = cp.sqrt(PixA_varREF_GPU + PixA_varSCI_GPU)
            PixA_DSNR_GPU = PixA_DCDIFF_GPU / PixA_NDIFF_GPU

        return PixA_DIFF_GPU, PixA_DCDIFF_GPU, PixA_DSNR_GPU
