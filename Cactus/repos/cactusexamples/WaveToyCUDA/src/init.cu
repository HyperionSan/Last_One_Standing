// -*-C++-*-

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cuda_runtime.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>



namespace WaveToyCUDA {
  
  // Check a return value, and if there is an error, output a
  // human-readable error message
  void check_error(cudaError_t cerr, char const *msg = "", ...)
#ifdef __GNUC__
  __attribute__((format (printf, 2, 3)))
#endif
  ;
  void check_error(cudaError_t cerr, char const *msg, ...)
  {
    if (cerr) {
      if (strcmp(msg, "")) {
        va_list ap;
        va_start(ap, msg);
        char *usermsg;
        vasprintf(&usermsg, msg, ap);
        va_end(ap);
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "CUDA Error %d: %s:\n%s",
                   (int)cerr, cudaGetErrorString(cerr), usermsg);
        free(usermsg);
      } else {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "CUDA Error %d: %s\n",
                   (int)cerr, cudaGetErrorString(cerr));
      }
    }
  }
  
  
  
  // Access a grid function in a kernel
  __device__ CCTK_REAL& gfelt(cudaPitchedPtr const& u,
                              size_t i, size_t j, size_t k)
  {
    return *(CCTK_REAL*)&((char*)u.ptr)[i*sizeof(CCTK_REAL) +
                                        j*u.pitch +
                                        k*u.pitch*u.ysize];
  }
  
  // Data living in the device memory
  namespace dev {
    cudaExtent ext;
    cudaPitchedPtr u;
  } // namespace dev
  
  
  
  // A simple kernel
  __global__ void init(size_t const lsh0, size_t const lsh1, size_t const lsh2,
                       cudaPitchedPtr const u)
  {
    size_t const i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t const j = blockIdx.y * blockDim.y + threadIdx.y;
    size_t const k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i<lsh0 && j<lsh1 && k<lsh2) {
      gfelt(u,i,j,k) = 0.0;
    }
  }
  
  
  
  extern "C"
  void WaveToyCUDA_Init(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    cudaError_t cerr;
    
    
    
    // Output device properties
    cudaDeviceProp prop;
    cerr = cudaGetDeviceProperties(&prop, 0);
    check_error(cerr, "Could not get device properties");
    printf("CUDA device properties (device %d):\n", 0);
    printf("   name:                        %s\n", prop.name);
    printf("   totalGlobalMem:              %zu\n", prop.totalGlobalMem);
    printf("   sharedMemPerBlock:           %zu\n", prop.sharedMemPerBlock);
    printf("   regsPerBlock:                %d\n", prop.regsPerBlock);
    printf("   warpSize:                    %d\n", prop.warpSize);
    printf("   memPitch:                    %zu\n", prop.memPitch);
    printf("   maxThreadsPerBlock:          %d\n", prop. maxThreadsPerBlock);
    printf("   maxThreadsDim:               %d %d %d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    printf("   maxGridSize:                 %d %d %d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    printf("   clockRate:                   %d\n", prop.clockRate);
    printf("   totalConstMem                %zu\n", prop.totalConstMem);
    printf("   major:                       %d\n", prop.major);
    printf("   minor:                       %d\n", prop.minor);
    printf("   textureAlignment:            %zu\n", prop.textureAlignment);
    printf("   deviceOverlap:               %d\n", prop.deviceOverlap);
    printf("   multiProcessorCount:         %d\n", prop.multiProcessorCount);
    printf("   kernelExecTimeoutEnabled:    %d\n", prop.kernelExecTimeoutEnabled);
    printf("   integrated:                  %d\n", prop.integrated);
    printf("   canMapHostMemory:            %d\n", prop.canMapHostMemory);
    printf("   computeMode:                 %d\n", prop.computeMode);
    printf("   maxTexture1D:                %d\n", prop.maxTexture1D);
    printf("   maxTexture2D:                %d %d\n", prop.maxTexture2D[0], prop.maxTexture2D[1]);
    printf("   maxTexture3D:                %d %d %d\n", prop.maxTexture3D[0], prop.maxTexture3D[1], prop.maxTexture3D[2]);
    printf("   maxTexture1DLayered:         %d %d\n", prop.maxTexture1DLayered[0], prop.maxTexture1DLayered[1]);
    printf("   maxTexture2DLayered:         %d %d %d\n", prop.maxTexture2DLayered[0], prop.maxTexture2DLayered[1], prop.maxTexture2DLayered[2]);
    printf("   surfaceAlignment:            %zu\n", prop.surfaceAlignment);
    printf("   concurrentKernels:           %d\n", prop.concurrentKernels);
    printf("   ECCEnabled:                  %d\n", prop.ECCEnabled);
    printf("   pciBusID:                    %d\n", prop.pciBusID);
    printf("   pciDeviceID:                 %d\n", prop.pciDeviceID);
    printf("   pciDomainID:                 %d\n", prop.pciDomainID);
    printf("   tccDriver:                   %d\n", prop.tccDriver);
    printf("   asyncEngineCount:            %d\n", prop.asyncEngineCount);
    printf("   unifiedAddressing:           %d\n", prop.unifiedAddressing);
    printf("   memoryClockRate:             %d\n", prop.memoryClockRate);
    printf("   memoryBusWidth:              %d\n", prop.memoryBusWidth);
    printf("   l2CacheSize:                 %d\n", prop.l2CacheSize);
    printf("   maxThreadsPerMultiProcessor: %d\n", prop.maxThreadsPerMultiProcessor);
    
    
    
    // Allocate memory
    if (!dev::u.ptr) {
      dev::ext = make_cudaExtent
        (sizeof(CCTK_REAL)*cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
      cerr = cudaMalloc3D(&dev::u, dev::ext);
      check_error(cerr, "Failed to allocate [%d,%d,%d] array",
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
    }
    
    // Call a kernel
    dim3 const blockDim(4, 4, 4);
    dim3 const gridDim((cctk_lsh[0] + blockDim.x - 1) / blockDim.x,
                       (cctk_lsh[1] + blockDim.y - 1) / blockDim.y,
                       (cctk_lsh[2] + blockDim.z - 1) / blockDim.z);
    init<<<gridDim, blockDim>>>(cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
                                dev::u);
    cerr = cudaGetLastError();
    check_error(cerr, "Could not call routine \"init\"");
    
    // Copy data to host
    cudaMemcpy3DParms parms = {0};
    parms.srcPtr = dev::u;
    parms.dstPtr = make_cudaPitchedPtr
      (u, sizeof(CCTK_REAL)*cctk_ash[0], cctk_lsh[0], cctk_lsh[1]);
    parms.extent = dev::ext;
    parms.kind = cudaMemcpyDeviceToHost;
    cerr = cudaMemcpy3D(&parms);
    check_error(cerr, "Failed to copy [%d,%d,%d] array",
                cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
    
    // Output grid function
    // for (ptrdiff_t k=0; k<cctk_lsh[2]; ++k) {
    //   for (ptrdiff_t j=0; j<cctk_lsh[1]; ++j) {
    //     for (ptrdiff_t i=0; i<cctk_lsh[0]; ++i) {
    //       ptrdiff_t const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
    //       printf("[%td,%td,%td]=%g\n", i,j,k, (double)u[ind3d]);
    //     }
    //   }
    // }
  }
  
} // namespace WaveToyCUDA
