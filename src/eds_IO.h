////////////////////////////////////////////////////////////////////////////////
//
// eds_IO.h
//
// created November 29, 2010 by HJvE
// Last Modified May 19, 2011 by HJvE
//
////////////////////////////////////////////////////////////////////////////////
#ifndef EDS_IO_H_
#define EDS_IO_H_

#include "environment.h"
#include <string.h>

#include "eds_2D_regular.h"
#include "extramath.h"
#include "hdf5.h"

////////////////////////////////////////////////////////////////////////////////

class c_eds_IO : public virtual c_eds
{
  public:
  
  bool save_EDS;
  
  // additional routines
  c_eds_IO();
  ~c_eds_IO();
  void save(const char *filenamebase, int entry);
};  

#endif // EDS_IO_2D_H_
