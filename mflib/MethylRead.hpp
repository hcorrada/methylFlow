#include <vector>
#include <string>
#include <lemon/list_graph.h>

#ifndef METHYLREAD_H
#define METHYLREAD_H

namespace methylFlow {

  enum ReadComparison { IDENTICAL, SUBREAD, SUPERREAD, METHOVERLAP, OVERLAP, NONE };

class MethylRead {
public:
  friend class MFGraph;

  MethylRead(int start, int length);
  MethylRead(const MethylRead &read);
  ~MethylRead();
  
  bool isMethConsistent(MethylRead *other);
  ReadComparison compare(MethylRead *other);
  int parseMethyl(std::string methylString);
  int merge(MethylRead *other);

  const std::string getMethString() const;
  std::string getMethString();

  std::string getString();
  const int start() const;
  const int length() const;
  const int end() const;
  const std::size_t ncpgs() const;

protected:
  // TODO: we need to distinguish region coordinates for modeling and read coordinates for genome coverage
  int rPos, rLen;
  std::vector<int> cpgOffset;
  std::vector<bool> methyl;
  int coverage;
  lemon::ListDigraph::Node node;
};

  inline const int MethylRead::start() const
  {
    return rPos;
  }

  inline const int MethylRead::length() const
  {
    return rLen;
  }

  // this assumes positions are 1-index
  inline const int MethylRead::end() const
  {
    return rPos + rLen - 1;
  }

  inline const std::size_t MethylRead::ncpgs() const
  {
    return this->cpgOffset.size();
  }

  struct CompareReadStarts : public std::binary_function<MethylRead *, MethylRead *, bool>
  {
  public:
    bool operator()(const MethylRead *x, const MethylRead *y) const;
  };

} // namespace methylFlow

#endif // METHYLREAD_H
