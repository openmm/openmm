#ifdef PLATFORM_Linux
#if defined(__i386__) || defined(__x86_64__)
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
#endif
#endif
#ifndef XTC
#define XTC 1
#include"xdrfile.h"
typedef struct float4
{
  float x, y, z, w;
} float4;

typedef struct double4
{
  double x, y, z, w;
} double4;

#ifdef __cplusplus
extern "C"
{
#endif

  // struct XTC_frame *xtc_read_frame(char *filename, int *natoms, int frame);
  void xtc_read_frame(char *filename, float *coords_arr, float *box_arr, float *time_arr, int *step_arr, int natoms, int frame, int nframes, int fidx);

  int xtc_write(char *filename, int natoms, int nframes, int *step, float *timex, float *pos, float *box);

  int xtc_write_frame(char *filename, int natoms, int nframes, int *step, float *timex, float *pos, float *box);
  // int xtc_write( char *filename, int natoms, int step, float time, float *pos, float*  box ) ;

  int xtc_truncate_to_step(char *infile, unsigned long maxstep);

  struct XTC_frame
  {
    float box[3];
    int natoms;
    unsigned long step;
    double time;
    float *pos;
  }; // XTC_frame;

  struct XTC_frame *xtc_read(char *filename, int *natoms, int *Nframes, double *dt, int *dstep);

  void xtc_read_new(char *filename, float *coords_arr, float *box_arr, float *time_arr, int *step_arr, int natoms, int nframes);

  int xtc_nframes(char *filename);
  int xtc_natoms(char *filename);


  XDRFILE* xtc_open_for_write(char* filename);
  int xtc_write_frames(XDRFILE* file, int natoms, int nframes, int *step, float *timex, float *pos, float *box);

#ifdef __cplusplus
}
#endif

#endif
