////////////////////////////////////////////////////////////////////////////////
//
// param_IO.cpp
//
// Created pre July 28, 2010 by HJvE
// Last modified September 7, 2011 by HJvE
//
// This file contains all I/O routines dealing with text parameter files.
//
////////////////////////////////////////////////////////////////////////////////

#include "param_IO.h"

////////////////////////////////////////////////////////////////////////////////

bool file_exists(const char* filename)
{
  FILE *pfile;
  pfile = fopen(filename,"r");
  if (pfile == NULL) return false;
  fclose(pfile);
  return true;
}

///////////////////////////////////////////////////////////////////////////////
// I/O functions for parameter files (parfiles)
////////////////////////////////////////////////////////////////////////////////

bool seek_parfile_par(FILE * pfile, const char* parameter)
{
  int x, xold = '\n';

  // printf("GOING TO SEEK THE FOLLOWING PARAMETER: %s\n", parameter);

  int parlength = strlen(parameter); // string length
  int i=0; // number of matching characters

  bool parfound = false;
  rewind(pfile);

  while(!parfound && !feof(pfile))
  {
    x = fgetc(pfile); 
    switch(x)
    {
      case '!': 
        while(x!='\n' && !feof(pfile)) 
        { 
          x=fgetc(pfile); 
          // printf("%c", x); fflush(stdout);
        }
        i=0; 
      break;
      default:
        i++;
        if (i > parlength)
        {
          if (x==' ' or x=='=') { parfound=true; }
            else i = 0;
        }
        else
          if (i == 1)
          {
            if (x!=parameter[i-1] or (xold != ',' and xold != '\n' and
              xold != '.' and xold != ';' and xold != 9 and xold != ' ')) i = 0;
          }
          else
          {
            if (x!=parameter[i-1]) i = 0;
          }
      break;  
    }
    xold = x;
  }
  if (parfound) return true; else return false;
}

////////////////////////////////////////////////////////////////////////////////

double get_parfile_double(FILE * pfile)
{
  size_t size_check;
  int c;
  double x = 0.0;

  // forward to first number
  do { c = fgetc(pfile); } while (!(c > 47 && c< 58) && c!='.' && c!='-');
  
  // one step backwards to reread the first number
  fseek(pfile, -1, SEEK_CUR);
  size_check = fscanf(pfile,"%lf",&x);
  
  if (size_check != 1)
  {
    printf("param_IO. get_parfile_double : read error\n");
    fflush(stdout);
    abort();
  }

  return x;
}

////////////////////////////////////////////////////////////////////////////////

int get_parfile_int(FILE * pfile)
{
  int c;
  int x = 0;
  size_t size_check;

  // forward to first number
  do {c = fgetc(pfile);} while (!(c > 47 && c< 58) && c!='.' && c!='-');

  // one step backwards to reread the first number
  fseek(pfile, -1, SEEK_CUR);
  size_check = fscanf(pfile,"%d",&x);
  
  if (size_check != 1)
  {
    printf("param_IO. get_parfile_double : read error\n");
    fflush(stdout);
    abort();
  }

  return x;
}

////////////////////////////////////////////////////////////////////////////////

long get_parfile_long(FILE * pfile)
{
  size_t size_check;
  int c;
  long x = 0;
  
  // forward to first number
  do {c = fgetc(pfile);} while (!(c > 47 && c< 58) && c!='.' && c!='-');

  // one step backwards to reread the first number
  fseek(pfile, -1, SEEK_CUR);
  size_check = fscanf(pfile,"%ld",&x);
  
  if (size_check != 1)
  {
    printf("param_IO. get_parfile_double : read error\n");
    fflush(stdout);
    abort();
  }
  return x;
}

////////////////////////////////////////////////////////////////////////////////

void get_parfile_string(FILE * pfile, char * &x)
{
  int c;
  int stringlength=0;
  int i;
  char* x1;
  char *x2; // strings WITHOUT \0 character
  bool finished=false;

  // seek beginning of string
  do 
  {
    c = fgetc(pfile);
  } 
  while (c == ' ' or c==9 or c=='=' or c=='"' or c=='\'');
  
  // save first character
  x1 = new char[1]; 
  x1[0] = (char) c;
  stringlength = 1;

  // get the whole string
  do
  {
    c=fgetc(pfile);
    if (c=='\n' or feof(pfile) or c==' ' or c==9 or c==',' or c==';'
      or c=='"' or c=='\'') finished = true;
    if (!finished)
    {
      stringlength++;
      x2 = new char[stringlength];
      for (i = 0; i < stringlength - 1; i++) x2[i] = x1[i];
      delete[] x1;
      x2[stringlength - 1] = (char) c;
      x1 = new char[stringlength];
      for (i = 0; i < stringlength; i++) x1[i] = x2[i];
      delete[] x2;
    }
  } while (!finished);
  delete[] x;

  x = new char[stringlength + 1];
  for (i = 0;i < stringlength; i++) x[i] = x1[i];
  x[stringlength] = '\0';
  delete[] x1;
}

////////////////////////////////////////////////////////////////////////////////

double double_from_parfile(char *filename, const char* parameter)
{
  double x;
  
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    abort();
  }

  FILE *pfile;
  pfile = fopen(filename,"r");
  if (!seek_parfile_par(pfile, parameter))
  {
    printf("Variable %s not found in %s! Using the value 1.0.\n",
      parameter, filename);
    fclose(pfile);
    x = 1.0;
  }
  else x = get_parfile_double(pfile);
  fclose (pfile);

  return x;
}

////////////////////////////////////////////////////////////////////////////////

int int_from_parfile(char *filename, const char* parameter)
{
  //printf("requesting int labeled %s from parfile %s\n", parameter, filename);
  
  int x;

  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    abort();
  }
  
  FILE *pfile;
  pfile = fopen(filename,"r");
  if (!seek_parfile_par(pfile, parameter))
  {
    printf("Variable %s not found in %s! Using the value 1.\n",
      parameter, filename);
    fclose(pfile);
    x = 1;
  }
  else x = get_parfile_int(pfile);
  fclose (pfile);

  return x;
}

////////////////////////////////////////////////////////////////////////////////

long long_from_parfile(char *filename, const char* parameter)
{
  long x;

  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    abort();
  }

  FILE *pfile;
  pfile = fopen(filename,"r");
  if (!seek_parfile_par(pfile, parameter))
  {
    printf("Variable %s not found in %s! Using the value 1.\n",
      parameter, filename);
    fclose(pfile);
    x = 1;
  }
  else x = get_parfile_int(pfile);
  fclose (pfile);

  return x;
}

////////////////////////////////////////////////////////////////////////////////

void string_from_parfile(char * filename, const char* parameter, char * &x)
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    fflush(stdout);
    exit(1);
  }

  FILE *pfile;
  pfile = fopen(filename,"r");
  if (!seek_parfile_par(pfile, parameter))
  {
    printf("Variable %s not found in %s! String not modified.\n", parameter, 
      filename);
    fflush(stdout);
  }
  else 
  {
    get_parfile_string(pfile, x);
  }
  fclose(pfile);
}

////////////////////////////////////////////////////////////////////////////////

void double_series_from_parfile(char * filename, const char * parameter, int no, 
  double x[])
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    abort();
  }

  int i;
  FILE * file_p;
  file_p = fopen(filename,"r");
  if (!seek_parfile_par(file_p, parameter))
    printf("Variable %s not found in %s! No doubles written to array.\n", 
      parameter, filename);
  else
    for(i=0;i<no;i++) x[i] = get_parfile_double(file_p);
  fclose(file_p);
}

////////////////////////////////////////////////////////////////////////////////

void int_series_from_parfile(char * filename, const char * parameter, int no, 
  int x[])
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    abort();
  }

  int i;
  FILE * file_p;
  file_p = fopen(filename,"r");
  if (!seek_parfile_par(file_p, parameter))
    printf("Variable %s not found in %s! No integers written to array.\n", 
      parameter, filename);
  else
    for(i=0;i<no;i++) x[i] = get_parfile_int(file_p);

  fclose(file_p);
}

////////////////////////////////////////////////////////////////////////////////

bool read_parfile_string(char *filename, const char* parameter, char * &x)
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    fflush(stdout);
    return false;
  }

  FILE *pfile;
  pfile = fopen(filename,"r");
  if (!seek_parfile_par(pfile, parameter))
  {
    fclose(pfile);
    return false;
  }
  else 
  {
    get_parfile_string(pfile, x);
    fclose(pfile);
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool read_parfile_double(char *filename, const char* parameter, double &x)
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    fflush(stdout);
    return false;
  }

  FILE * file_p;
  file_p = fopen(filename,"r");
  if (!seek_parfile_par(file_p, parameter))
  {
    x = -1.; // default value
    return false;
  }
  else
    x = get_parfile_double(file_p);

  fclose(file_p);
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool read_parfile_int(char *filename, const char* parameter, int &x)
{
  if (!file_exists(filename))
  {
    printf("ERROR: file %s does not exist.\n", filename);
    fflush(stdout);
    return false;
  }

  FILE * file_p;
  file_p = fopen(filename,"r");
  if (!seek_parfile_par(file_p, parameter))
  {
    x = -1; // default value
    return false;
  }
  else
    x = get_parfile_int(file_p);

  fclose(file_p);
  
  return true;
}
