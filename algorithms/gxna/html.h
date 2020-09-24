#ifndef __HTML_H
#define __HTML_H

// classes and functions to write HTML
// very basic, little error checking performed

#include <iostream>
using namespace std;

class htmlstream {
  ostream str;
 public:
  void title(const string&);
  void twoFrames(const string&, const string&, double);
  void beginTable(int border = 1);
  void endTable(void);
};

#endif
