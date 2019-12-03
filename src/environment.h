////////////////////////////////////////////////////////////////////////////////
//
// environment.h
//
// Created: July 28, 2010, by HJvE
// Last modified: June 8, 2016 by HJvE
//
// All pre-compiler settings are stored in this file.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

////////////////////////////////////////////////////////////////////////////////

// include the following variables in you main program file
extern int numprocs; // stores the number of processors used for this run
extern int myid; // identity of single processor

////////////////////////////////////////////////////////////////////////////////

#define ENABLED_ 1
#define DISABLED_ 2

// settings to implement. Change these for different implementations
#define OPEN_MPI_           ENABLED_
  // if OPEN_MPI_ is disabled, openMP will be used occasionally instead of MPI
  // and shared memory between cores is utilized. 
#define BOOST_              DISABLED_
  // enabled for a Lorentz-boosted frame BOX / simulation

//------------------------------------------------------------------------------
// double checks to prevent typo's above from leading to hard-to-diagnose issues

#if BOOST_ != ENABLED_
#if BOOST_ != DISABLED_
#error typo in BOOST_ setting. ENABLED_ or DISABLED_ spelled wrongly
#endif
#endif

#ifndef BOOST_
#error Typo in BOOST_
#endif

#if OPEN_MPI_ != ENABLED_
#if OPEN_MPI_ != DISABLED_
#error typo in OPEN_MPI_ setting. ENABLED_ or DISABLED_ spelled wrongly
#endif
#endif

#ifndef OPEN_MPI_
#error Typo in OPEN_MPI_
#endif

///////////////////////////////////////////////////////////////////////////////

// Leave at DISABLED for reliable early-time self-absorption results, especially
// in the wind case
#define LOCAL_SELF_ABSORPTION_  DISABLED_

#endif // ENVIRONMENT_H_
