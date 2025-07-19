Functional Connectivity Gradient Utilities

Core functions

FCG_OrthProc.m - Orthogonal Procrustes registration: estimate rigid rotation, one source gradient to a target
FCG_OrthProcReg.m - Same as FCG_OrthProcReg.m, but for set of source gradients to a target
FCG_Reg1.m - Industry standard - Orthogonal Procrustes of set to first, compute mean, register to mean
FCG_ApplyReg.m - Apply array of registrations (rigid plus optional scale) to set of gradients
FCG_AtlasMean.m - Apply array of registrations (rigid plus optional scale) to set of gradients, compute mean
FCG_ScaleReg.m - Estimate scale, one source gradient to a target
FCG_RegAll.m - For set of gradients, estimate all possible registrations, pick most typical as target;
               allows for orthogonal, scaled orthogonal or full affine

Utility functions

FCG_Load.m - Load gradient data
FCG_PlotMatrix.m - Specialiesd matrix plot function for gradients
FCG_Print.m - Full page print, adding page-wise title
FCG_PrintGrads.m - For a pair of gradients, multipage printing of every subject


Depricated

FCG_Reg.m - Combines FCG_Load, outlier filtering, and FCG_RegAll, and different similarity metrics
            to define the 'most typical' subject


