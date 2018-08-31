////////////////////////////////////////////////////////////////////////////////
//
// Param_IO.h
//
// Created pre July 28, 2010 by HJvE
// Last Modified September 8, 2011 by HJvE
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PARAM_IO_H_
#define PARAM_IO_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// General functions
////////////////////////////////////////////////////////////////////////////////

bool file_exists(const char* filename);
// returns true if file exist, false otherwise

////////////////////////////////////////////////////////////////////////////////
// I/O functions for parameter file interactions
////////////////////////////////////////////////////////////////////////////////

bool seek_parfile_par(FILE * pfile, const char* parameter);
// This function looks for a string parameter in the file pfile. It returns true
// if this string is present, and not part of a larger string (e.g. it checks
// the surrounding positions for spaces, comma's etc to distinguish for example
// between 'rho' and 'rho_1'). It also sets the file pointer to the 2nd position
// after the string / parameter name. It returns false if the parameter is not
// found.

double get_parfile_double(FILE * pfile);
// assuming the file pointer is at a position where the next following number
// is a double, read and return this double.

int get_parfile_int(FILE * pfile);

long get_parfile_long(FILE *pfile);

void get_parfile_string(FILE * pfile, char *&x);
// IMPORTANT: though the array x will be overwritten, it has to be defined
// and allocated using 'new'!

double double_from_parfile(char *filename, const char* parameter);
// In the file specified with filename, looks for a parameter specified with
// parameter and returns the number following that parameter in the file. It
// uses seek_parfile_par and get_parfile_double and opens and closes the file.
 
int int_from_parfile(char * filename, const char* parameter);

long long_from_parfile(char *filename, const char *parameter);

void string_from_parfile(char * filename, const char* parameter, char * &x);

void double_series_from_parfile(char * filename, const char * parameter, int no, 
  double x[]);

void int_series_from_parfile(char * filename, const char * parameter, int no, 
  int x[]);

////////////////////////////////////////////////////////////////////////////////
// Upgraded IO functions included boolean return value

bool read_parfile_string(char *filename, const char* parameter, char * &x);
// in the file specified with filename, look for a parameter specified with the
// string parameter and return the string following that parameter in the file
// into the variable x. The function returns true upon success

bool read_parfile_double(char *filename, const char* parameter, double &x);
// in the file specified with filename, look for a parameter specified with the
// string parameter and return the number following that parameter in the file
// into the variable x. The function returns true upon success

bool read_parfile_int(char *filename, const char* parameter, int &x);
// in the file specified with filename, look for a parameter specified with the
// string parameter and return the number following that parameter in the file
// into the variable x. The function returns true upon success

#endif /*PARAM_IO_H_*/

