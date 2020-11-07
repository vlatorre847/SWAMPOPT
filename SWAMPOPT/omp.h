#ifndef _OMP_H_DEF
#define _OMP_H_DEF

#include <stdlib.h> // support size_t

/*
   7 Updated to support OpenMP 4.5, B.1 page 327
   8
   9 Liao 10/10/2016
   10
   11 From the specification 3.0, Chapter 3.1, page 108
   12 What to put into this file:
   13 * The prototypes of all the runtime library routines, with "C" linkage
   14 * The type omp_lock_t
   15 * The type omp_nest_lock_t
   16 * The type omp_sched_t
   17
   18 Specification 3.0 page 302
   19 D.1 Example of the omp.h Header File
   20  * */
// This is not correct. omp.h is not required to trigger _OPENMP macro
//#define _OPENMP 201511
typedef void *omp_lock_t;  /* represented as a pointer */
typedef void *omp_nest_lock_t; /* represented as a pointer */

/*
    28 * define the lock hints
    29 */
typedef enum omp_lock_hint_t
{
       omp_lock_hint_none = 0,
       omp_lock_hint_uncontended = 1,
       omp_lock_hint_contended = 2,
       omp_lock_hint_nonspeculative = 4,
       omp_lock_hint_speculative = 8
         /* , Add vendor specific constants for lock hints here,
            38        starting from the most-significant bit. */
     } omp_lock_hint_t;


 /*
    43 *  define the schedule kinds
    44 */

 typedef enum omp_sched_t
 {
       omp_sched_static = 1,
       omp_sched_dynamic = 2,
       omp_sched_guided = 3,
       omp_sched_auto = 4
     }
omp_sched_t;

 /*
    55 * define the proc bind values
    56 */
 typedef enum omp_proc_bind_t
 {
       omp_proc_bind_false = 0,
       omp_proc_bind_true = 1,
       omp_proc_bind_master = 2,
       omp_proc_bind_close = 3,
       omp_proc_bind_spread = 4
     } omp_proc_bind_t;


 #ifdef __cplusplus
 extern "C" {
     #endif
    
     /*
        72  * Execution Environment Functions
        73  */
     extern void omp_set_num_threads(int num);
     extern int omp_get_num_threads(void);
     extern int omp_get_max_threads(void);
     extern int omp_get_thread_num(void);
     extern int omp_get_num_procs(void);
    
     int omp_in_parallel(void);
     void omp_set_dynamic(int dynamic_thds);
     int omp_get_dynamic(void);
    
     int omp_get_cancellation(void);
    
     void omp_set_nested(int n_nested);
     int omp_get_nested(void);
    
     /*
        90  * Other internal variables
        91  */
    
     void omp_set_schedule (omp_sched_t, int);
     void omp_get_schedule (omp_sched_t *, int *);
     int omp_get_thread_limit (void);
     void omp_set_max_active_levels (int);
     int omp_get_max_active_levels (void);
     int omp_get_level (void);
     int omp_get_ancestor_thread_num (int);
     int omp_get_team_size (int);
     int omp_get_active_level (void);
    
     int omp_in_final(void);
     omp_proc_bind_t omp_get_proc_bind(void);
     int omp_get_num_places(void);
     int omp_get_place_num_procs(int place_num);
     void omp_get_place_proc_ids(int place_num, int *ids);
     int omp_get_place_num(void);
     int omp_get_partition_num_places(void);
     void omp_get_partition_place_nums(int *place_nums);
    
     /*
         113  * Support accelerators as target devices
         114  */
    
     void omp_set_default_device(int device_num);
     int omp_get_default_device(void);
    
     /*  find the max number of devices on the system */
     int omp_get_max_devices(void);
     /* set number of active devices to be used */
     void omp_set_num_devices(int);
    
     /* get number of available devices */
     int omp_get_num_devices(void);
     // GCC 4.0 provides omp_get_num_devices() already, but without supporting GPUs
     // I have to use another function to bypass it
     int xomp_get_num_devices(void);
    
     int omp_get_num_teams(void);
     int omp_get_team_num(void);
    
     int omp_is_initial_device(void);
     int omp_get_initial_device(void);
     int omp_get_max_task_priority(void);
    
     /*
         138  * Lock Functions
         139  */
     void omp_init_lock(omp_lock_t *lock);
     void omp_init_lock_with_hint(omp_lock_t *lock,
                                                                          omp_lock_hint_t hint);
     void omp_destroy_lock(omp_lock_t *lock);
     void omp_set_lock(omp_lock_t *lock);
     void omp_unset_lock(omp_lock_t *lock);
     int omp_test_lock(omp_lock_t *lock);
     void omp_init_nest_lock(omp_nest_lock_t *lock);
     void omp_init_nest_lock_with_hint(omp_nest_lock_t *lock,
                                                                                    omp_lock_hint_t hint);
     void omp_destroy_nest_lock(omp_nest_lock_t *lock);
     void omp_set_nest_lock(omp_nest_lock_t *lock);
     void omp_unset_nest_lock(omp_nest_lock_t *lock);
     int omp_test_nest_lock(omp_nest_lock_t *lock);
    
     /*
         156  * Timer routine
         157  */
     double omp_get_wtime(void);
     double omp_get_wtick(void);
    
     void * omp_target_alloc(size_t _size, int _device_num);
     void omp_target_free(void * device_ptr, int _device_num);
     int omp_target_is_present(void * ptr, int _device_num);
     int omp_target_memcpy(void *dst, void *src, size_t _length,
                                                     size_t dst_offset, size_t src_offset,
                                                     int dst_device_num, int src_device_num);
    
     int omp_target_memcpy_rect(
                                                               void *dst, void *src,
                                                               size_t element_size,
                                                               int num_dims,
                                                               const size_t *volume,
                                                               const size_t *dst_offsets,
                                                               const size_t *src_offsets,
                                                               const size_t *dst_dimensions,
                                                               const size_t *src_dimensions,
                                                               int dst_device_num, int src_device_num);
    
     int omp_target_associate_ptr(void * host_ptr,
                                                                   void * device_ptr,
                                                                   size_t _size,
                                                                   size_t device_offset,
                                                                   int device_num);
    
     int omp_target_disassociate_ptr(void * ptr,
                                                                         int device_num);
    
    
     /*
         190  * FORTRAN Execution Environment Function Wrappers
         191  * Fortran stuff should be handled in omp_lib.h
         192 void omp_set_num_threads_(int *num);
         193 int omp_get_num_threads_(void);
         194 int omp_get_max_threads_(void);
         195 int omp_get_thread_num_(void);
         196 int omp_get_num_procs_(void);
         197 int omp_in_parallel_(void);
         198 void omp_set_dynamic_(int *dynamic_thds);
         199 int omp_get_dynamic_(void);
         200 void omp_set_nested_(int *n_nested);
         201 int omp_get_nested_(void);
         202
         203  * FORTRAN Lock Function Wrappers
         204 typedef unsigned int _omf77Lock_t;
         205 void omp_init_lock_(_omf77Lock_t *lock);
         206 void omp_init_nest_lock_(_omf77Lock_t *lock);
         207 void omp_destroy_lock_(_omf77Lock_t *lock);
         208 void omp_destroy_nest_lock_(_omf77Lock_t *lock);
         209 void omp_set_lock_(_omf77Lock_t *lock);
         210 void omp_set_nest_lock_(_omf77Lock_t *lock);
         211 void omp_unset_lock_(_omf77Lock_t *lock);
         212 void omp_unset_nest_lock_(_omf77Lock_t *lock);
         213 int omp_test_lock_(_omf77Lock_t *lock);
         214 int omp_test_nest_lock_(_omf77Lock_t *lock);
         215
         216 */
     #ifdef __cplusplus
     } /* closing brace for extern "C" */
 #endif

 #endif /* _OMP_H_DEF */
