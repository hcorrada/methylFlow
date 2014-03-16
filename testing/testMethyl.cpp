#include "mflib/MethylRead.hpp"
#include <cassert>
#include <iostream>

int main() {
  methylFlow::MethylRead m1(3, 10);
  m1.parseMethyl("6:M,8:M");

  methylFlow::MethylRead m2(3, 10);
  m2.parseMethyl("6:M,8:M");
  assert(m1.compare(&m2) == methylFlow::IDENTICAL);

  methylFlow::MethylRead m3(3, 15);
  m3.parseMethyl("6:M,8:M");
  assert(m1.compare(&m3) == methylFlow::SUPERREAD);

  methylFlow::MethylRead m4(3, 8);
  m4.parseMethyl("6:M,8:M");
  assert(m1.compare(&m4) == methylFlow::SUBREAD);

  methylFlow::MethylRead m5(5, 10);
  m5.parseMethyl("4:M,6:M");
  assert(m1.compare(&m5) == methylFlow::METHOVERLAP);

  
  methylFlow::MethylRead m6(5, 10);
  m6.parseMethyl("4:U,6:U");
  assert(m1.compare(&m6) == methylFlow::OVERLAP);

  methylFlow::MethylRead m7(5, 10);
  m7.parseMethyl("4:M,6:U");
  assert(m1.compare(&m7) == methylFlow::OVERLAP);

  methylFlow::MethylRead m8(15, 10);
  m8.parseMethyl("4:M");
  assert(m1.compare(&m8) == methylFlow::NONE);


  methylFlow::MethylRead *u = new methylFlow::MethylRead(1, 10);
  u->parseMethyl("2:M,6:M");

  methylFlow::MethylRead *v = new methylFlow::MethylRead(3,10);
  v->parseMethyl("4:M,9:M");

  u->merge(v);
  std::cout << u->getString() << std::endl;
  
  delete u;
  delete v;

  u = new methylFlow::MethylRead(24, 10);
  u->parseMethyl("2:U,9:M");
  
  v = new methylFlow::MethylRead(26, 10);
  v->parseMethyl("7:M,10:M");
  assert(u->compare(v) == methylFlow::METHOVERLAP);

  return 0;
}
