#include <string>
#include <iostream>
#include <sstream>
#include <cassert>

#include "MethylRead.hpp"

namespace methylFlow {
  MethylRead::MethylRead(int pos, int len) : rPos(pos), rLen(len), cpgOffset(), methyl(), coverage(0)
  {
  }

  MethylRead::MethylRead(const MethylRead &read) : rPos(read.start()), rLen(read.length()), coverage(read.coverage)
  {
    this->cpgOffset = std::vector<int>(read.cpgOffset);
    this->methyl = std::vector<bool>(read.methyl);
  }

  MethylRead::~MethylRead()
  {
  }

  int MethylRead::parseMethyl(std::string methString)
  {
    std::size_t found;
    std::size_t curStringOffset = 0; 
    int curPos = 0;
    bool curMeth = false;
      if (methString.length() == 0) {
          std::cout << "Error parsing methylation string, There is no read" << std::endl;
          return -1;

      }
      char  lastChar = methString.at( methString.length() - 1 );
      if (lastChar != ',') {
          methString.append(1, ',');
      }
    while (curStringOffset < methString.length()) {
      found = methString.find(":", curStringOffset);

      if (found == std::string::npos) {
          std::cout << "Error parsing methylation string, : is not found" << std::endl;
	return -1;
      }

      curPos = atoi(methString.substr(curStringOffset, found - curStringOffset).c_str());

      curStringOffset = found + 1;
      found = methString.find(",", curStringOffset);

      if (found == std::string::npos) {
	std::cout << "Error parsing methylation string, , is not found" << std::endl;
	return -1;
      }
      curMeth = (methString.substr(curStringOffset, found - curStringOffset) == "M"); 
      cpgOffset.push_back(curPos);
      methyl.push_back(curMeth);
      curStringOffset = found + 1;
    }
    return 0;
  }
    
    
  int MethylRead::parseXMtag(std::string XM){
      std::size_t foundU, foundM, found;
      std::size_t curStringOffset = XM.find("XM:Z:", 0) + 5;
      
      int curPos = 0;
      bool curMeth = false;
      
      while (curStringOffset < XM.length()) {
          foundM = XM.find("Z", curStringOffset);
          foundU = XM.find("z", curStringOffset);
          
          found = std::min(foundM, foundU);
          
          if (found == std::string::npos) {
              break;
          }
          curPos = found - 5;
          curMeth = (XM[found] == 'Z');
          cpgOffset.push_back(curPos);
          methyl.push_back(curMeth);
          curStringOffset = found + 1;
      }
      return 0;
      
  }


  const std::string MethylRead::getMethString() const
  {
    if (ncpgs() == 0) {
      return "*";
    }

    std::stringstream out;
    for (std::size_t i = 0; i != cpgOffset.size(); ++i) {
      if (i != 0) {
	out << ",";
      }
      out << cpgOffset[i] << ":" << (methyl[i] ? "M" : "U");
    }
    return out.str();
  }

  std::string MethylRead::getMethString()
  {
    const MethylRead * obj = const_cast<MethylRead *>(this);
    const std::string tmp = obj->getMethString();
    const std::string *p1 = &tmp;
    std::string *p2 = const_cast<std::string *>(p1);
    return *p2;
  }

  std::string MethylRead::getString()
  {
    std::stringstream out;
    out << "start: " << this->start();
    out << " end: " << this->end();
    out << " meth: " << this->getMethString();
    return out.str();
  }

    float MethylRead::distance(MethylRead* other, int &common) {
        int offset = other->start() - this->start();
        std::size_t j = 0;
        float match = 0;
        for (std::size_t i=0; i < this->cpgOffset.size(); ++i) {
            // advance pointer of other as long as needed
            while( j < other->cpgOffset.size() && (this->cpgOffset[i] - offset) > other->cpgOffset[j] )
                j++;
            
         //   std::cout << i << " " << j << std::endl;
         //   std::cout << this->cpgOffset.size() << " " << other->cpgOffset.size() << std::endl;
            // no more cpgs on other, so return true
            if (j == this->cpgOffset.size()) break;
            
            // check if pointers at same position
            if ( this->cpgOffset[i] - offset == other->cpgOffset[j] ) {
                // we are, check if consistent
                common++;
                if ( this->methyl[i] == other->methyl[j] ) match++;
            }
            
            // no more cpgs on other, so return true
            if (j == this->cpgOffset.size()) break;
        }
        
        //std::cout << "match 1 = " << match << std::endl;
        //std::cout << "match 2 = " << this->cpgOffset.size() << std::endl;
        
        if (this->cpgOffset.size() > 0)
            return match;
        else
            return 0;
        
    }
    
    void MethylRead::write() {
        std::cout << "write Methyl" << std::endl;
        for (unsigned int i=0; i<this->cpgOffset.size(); i++) {
            if (!this->methyl[i]) {
                std::cout << this->cpgOffset.at(i) << ":" << "U, ";
            } else{
                std::cout << this->cpgOffset.at(i) << ":" << "M, ";
            }
        }
        std::cout << std::endl;
    }
    
  bool MethylRead::isMethConsistent(MethylRead *other)
  {
    // always assume comparing left read to right read
    if (this->start() > other->start())
      return false;

    if (this->ncpgs() == 0 || other->ncpgs() == 0)
      return true;

    int offset = other->start() - this->start();
    std::size_t j = 0;
    for (std::size_t i=0; i < this->cpgOffset.size(); ++i) {
      // advance pointer of other as long as needed
      while( j < other->cpgOffset.size() && (this->cpgOffset[i] - offset) > other->cpgOffset[j] )
	j++;

      // no more cpgs on other, so return true
      if (j == this->cpgOffset.size()) return true;

      // check if pointers at same position
      if ( this->cpgOffset[i] - offset == other->cpgOffset[j] ) {
	// we are, check if consistent
	if ( this->methyl[i] != other->methyl[j] ) return false;
      }
    }
    return true;
  }

  ReadComparison MethylRead::compare(MethylRead *other) {
    if (other->start() > this->end())
      return NONE;

    if (!isMethConsistent(other))
      return OVERLAP;

    if (other->start() > this->start())
      return other->end() <= this->end() ? SUBREAD : METHOVERLAP;

    assert(other->start() == this->start());

    if (other->end() == this->end())
      return IDENTICAL;

    if (other->end() > this->end())
      return SUPERREAD;

    return SUBREAD; 
  }

  int MethylRead::merge(MethylRead *other)
  {
    int offset = other->start() - this->start();

    std::size_t j = 0;
    for (std::vector<int>::iterator i = this->cpgOffset.begin(); i != this->cpgOffset.end() && j < other->cpgOffset.size(); ++i) {
      if (*i == (other->cpgOffset[j] + offset) ) {
	++j;
      }
    }

    for (; j < other->cpgOffset.size(); ++j) {
      this->cpgOffset.push_back( other->cpgOffset[j] + offset );
      this->methyl.push_back( other->methyl[j] );
    }
    this->rLen = offset + other->rLen;
    return 0;
  }

  bool CompareReadStarts::operator()(const MethylRead *x, const MethylRead *y) const
    {
      return x->start() < y->start();
    }
  

} // namespace methylFlow
