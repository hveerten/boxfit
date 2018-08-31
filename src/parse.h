////////////////////////////////////////////////////////////////////////////////
//
// parse.h
//
// Created August 3, 2010 by HJvE
// Last Modified Dec 15, 2014 by HJvE
//
// set of tools to parse command line arguments
//
////////////////////////////////////////////////////////////////////////////////
#ifndef PARSE_H_
#define PARSE_H_

bool parse_int(const char *argname, int &x, int argc, char* argv[]);
// parses the command line arguments (argc indicates how many of them there are, 
// while argv is an array of strings containing the arguments themselves), for
// any to start with the string given by argname. The integer number following
// that string is stored in x. The function returns true if argname is found
// and false if argname is not found.
// equal signs etc need to be included in argname. For example, suppose 
// some number needs to be indicated by "-n=5", then argname = "-n=". No spaces
// are allowed. (because, for example, "-n = 5" would be stored as 3 separate
// arguments)

bool parse_float(const char *argname, float &x, int argc, char* argv[]);

bool parse_double(const char *argname, double &x, int argc, char* argv[]);

bool parse_string(const char *argname, char *&s, int argc, char* argv[]);
// note that the string s needs to be declared in memory. It will be deleted
// and redeclared matching sufficient length. Note that argname should 
// explicitly include "=", for exampe "-resolution="

bool parse(const char *argname, int argc, char* argv[]);
// returns true if argname provided, false otherwise
#endif // PARSE_H_
