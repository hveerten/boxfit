////////////////////////////////////////////////////////////////////////////////
//
// parse.cpp
//
// Created August 3, 2010 by HJvE
// Last modified Dec 15, 2014 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#include "parse.h"
#include <stdio.h>
#include <string.h>

bool parse_int(const char *argname, int &x, int argc, char* argv[])
{
  int i;
  int argno = -1;
  
  for (i=1; i < argc; i++)
  {
    if (strncmp(argv[i], argname, strlen(argname)) == 0) argno = i;
  }
  
  if (argno == -1) { x = -1; return false; }
  
  sscanf(argv[argno]+strlen(argname), "%d", &x);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool parse_float(const char *argname, float &x, int argc, char* argv[])
{
  int i;
  int argno = -1;
  
  for (i=1; i < argc; i++)
  {
    if (strncmp(argv[i], argname, strlen(argname)) == 0) argno = i;
  }
  
  if (argno == -1) { x = 0.0; return false; }
  
  sscanf(argv[argno]+strlen(argname), "%e", &x);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool parse_double(const char *argname, double &x, int argc, char* argv[])
{
  int i;
  int argno = -1;
  
  for (i=1; i < argc; i++)
  {
    if (strncmp(argv[i], argname, strlen(argname)) == 0) argno = i;
  }
  
  if (argno == -1) { x = 0.0; return false; }
  
  sscanf(argv[argno]+strlen(argname), "%le", &x);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool parse_string(const char *argname, char *&s, int argc, char* argv[])
{
  int i;
  int argno = -1;
  
  for (i=1; i < argc; i++)
  {
    if (strncmp(argv[i], argname, strlen(argname)) == 0) argno = i;
  }
  
  if (argno == -1) { return false; }
  
  // delete dummy declaration and declare with exactly the right length
  delete[] s;
  s = new char[strlen(argv[argno]) - strlen(argname) + 1];

  sscanf(argv[argno] + strlen(argname), "%s", s);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool parse(const char *argname, int argc, char* argv[])
{
  int i;
  int argno = -1;
  
  for (i=1; i < argc; i++)
  {
    if (strncmp(argv[i], argname, strlen(argname)) == 0) argno = i;
  }
  
  if (argno == -1) { return false; }

  return true;  
}
