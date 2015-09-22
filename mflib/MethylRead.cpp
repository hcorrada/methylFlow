#include <string>
#include <iostream>
#include <sstream>
#include <cassert>

#include "MethylRead.hpp"

namespace methylFlow {
  MethylRead::MethylRead(int pos, int len) : rPos(pos), rLen(len), cpgs(), coverage(0)
  {
  }

  MethylRead::MethylRead(const MethylRead &read) : rPos(read.start()), rLen(read.length()), coverage(read.coverage)
  {
    this->cpgs = std::vector<CpgEntry>(read.cpgs);
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
      CpgEntry entry = { curPos, curMeth };
      cpgs.push_back(entry);
      curStringOffset = found + 1;
    }
    return 0;
  }
    
    int MethylRead::parseXMtag(std::string XM, std::string CIGAR){
        //std::cout << "XM Test CIGAR: " << CIGAR << std::endl;
        //std::cout << "XM Test before: " << XM << std::endl;
        std::string revisedXM = reviseXMtag(XM, CIGAR);
        //std::cout << "XM Test after : " << revisedXM << std::endl;
        return parseRevisedXMtag(revisedXM);
        
    }
    
    
    std::string MethylRead::reviseXMtag(std::string XM, std::string CIGAR){
        if (CIGAR == "*")
            return XM;
        
        std::vector<CigarEntry> cigarEntries;
        std::size_t curStringOffset = 0;
        //We set start=5 because there is "XM:Z:" at the beginning of XM-Tag ...xh.......Z
        int start = 5;
        std::size_t found, foundM, foundD, foundI;
        while (curStringOffset < CIGAR.length()) {
            foundM = CIGAR.find("M", curStringOffset);
            foundD = CIGAR.find("D", curStringOffset);
            foundI = CIGAR.find("I", curStringOffset);
            
            found = std::min(std::min(foundM, foundD), foundI);
            
            if (found == std::string::npos) {
                break;
            }
            int length = std::atoi(CIGAR.substr(curStringOffset, found - curStringOffset).c_str());
            
            CigarEntry entry = {start, length, CIGAR[found]};
            cigarEntries.push_back(entry);
            
            if (CIGAR[found] == 'M' || CIGAR[found] == 'I')
                start += length;
                                 
            curStringOffset = found + 1;
        }
        
        for (int i = cigarEntries.size()-1; i >= 0; i--) {
            if (cigarEntries.at(i).indel == 'D') {
                XM.insert(cigarEntries.at(i).start, cigarEntries.at(i).length, '.');
            } else if (cigarEntries.at(i).indel == 'I') {
                XM.erase(cigarEntries.at(i).start, cigarEntries.at(i).length);
            }
        }
        return XM;
    }
    
    int MethylRead::parseRevisedXMtag(std::string XM){
    
      std::size_t foundU, foundM, found;
//      std::size_t curStringOffset = XM.find("XM:Z:", 0) + 5;
      std::size_t curStringOffset = 0;
      XM.erase(0, 5);
      
      int curPos = 0;
      bool curMeth = false;
      
      while (curStringOffset < XM.length()) {
          foundM = XM.find("Z", curStringOffset);
          foundU = XM.find("z", curStringOffset);
          
          found = std::min(foundM, foundU);
          
          if (found == std::string::npos) {
              break;
          }
          curPos = found;
          curMeth = (XM[found] == 'Z');

	      CpgEntry entry = { curPos, curMeth };
          cpgs.push_back(entry);
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
    for (std::vector<CpgEntry>::const_iterator it = this->cpgs.begin(); it != this->cpgs.end(); ++it) {
      if (it != this->cpgs.begin()) {
	out << ",";
      }
      CpgEntry entry = *it;
      out << entry.offset << ":" << (entry.methyl ? "M" : "U");
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

	for (std::size_t i = 0; i < this->cpgs.size(); ++i) {
	  CpgEntry thisEntry = this->cpgs[i];

            // advance pointer of other as long as needed
	  while( j < other->cpgs.size() && (thisEntry.offset - offset) > other->cpgs[j].offset )
                ++j;
            
         //   std::cout << i << " " << j << std::endl;
         //   std::cout << this->cpgOffset.size() << " " << other->cpgOffset.size() << std::endl;
            // no more cpgs on other, so return true
            if (j == this->cpgs.size()) break;
            
            // check if pointers at same position
	    CpgEntry otherEntry = other->cpgs[j];
            if ( thisEntry.offset - offset == otherEntry.offset ) {
                // we are, check if consistent
                common++;
                if ( thisEntry.methyl == otherEntry.methyl ) match++;
            }            
        }
        
        //std::cout << "match 1 = " << match << std::endl;
        //std::cout << "match 2 = " << this->cpgOffset.size() << std::endl;
        
        if (cpgs.size() > 0)
            return match;
        else
            return 0;
        
    }
    
    void MethylRead::write() {
        std::cout << "write Methyl" << std::endl;
        std::cout << "start = " <<  this->start() << std::endl;

	for (std::vector<CpgEntry>::iterator it = cpgs.begin(); it != cpgs.end(); ++it) {
	  CpgEntry entry = *it;
	  std::cout << entry.offset << ":" << (entry.methyl ? "M," : "U,");
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

    for (std::size_t i=0; i < this->cpgs.size(); ++i) {
      CpgEntry thisEntry = this->cpgs[i];

      // advance pointer of other as long as needed
      while( j < other->cpgs.size() && (thisEntry.offset - offset) > other->cpgs[j].offset )
	j++;

      // no more cpgs on other, so return true
      if (j == this->cpgs.size()) return true;

      // check if pointers at same position
      CpgEntry otherEntry = other->cpgs[j];
      if ( thisEntry.offset - offset == otherEntry.offset ) {
	// we are, check if consistent
	if ( thisEntry.methyl != otherEntry.methyl ) return false;
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
    for (std::vector<CpgEntry>::iterator i = this->cpgs.begin(); i != this->cpgs.end() && j < other->cpgs.size(); ++i) {
      if (i->offset == (other->cpgs[j].offset + offset) ) {
	++j;
      }
    }

    for (; j < other->cpgs.size(); ++j) {
      CpgEntry newEntry = other->cpgs[j];
      newEntry.offset += offset;
      this->cpgs.push_back( newEntry );
    }
    this->rLen = offset + other->rLen;
    return 0;
  }

  bool CompareReadStarts::operator()(const MethylRead *x, const MethylRead *y) const
    {
      return x->start() < y->start();
    }
  

} // namespace methylFlow
