#include "html.h"

void htmlstream::title(const string& title) {
  str << "<TITLE>" << title << "</TITLE>\n";
}

void htmlstream::twoFrames(const string& file1, const string& file2, double ratio) {
  int p1 = int(100 * (ratio / (1 + ratio)));
  int p2 = 100 - p1; 
  str << "<FRAMESET COLS=\"" << p1 << "%," << p2 << "%\">\n";
  str << "<FRAME SRC=\"" << file1 << "\">\n";
  str << "<FRAME SRC=\"" << file2 << "\">\n";
  str << "</FRAMESET>\n";
}

void htmlstream::beginTable(int border) {
  str << "<TABLE BORDER=\"" << border << "\">\n";
}

void htmlstream::endTable(void) {
  str << "</TABLE\n";
}
