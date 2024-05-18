#ifndef COOT_HEADER_INFO_HH
#define COOT_HEADER_INFO_HH

#include <string>
#include <vector>

namespace moorhen {

   class helix_t {
   public:
      int serNum;
      std::string helixID;
      std::string initResName; // name of the helix's initial residue
      std::string initChainID; // chain ID for the chain containing the helix
      int     initSeqNum;      // sequence number of the initial residue
      std::string initICode;   // insertion code of the initial residue
      std::string endResName;  // name of the helix's terminal residue
      std::string endChainID;  // chain ID for the chain containing the helix
      int     endSeqNum;       // sequence number of the terminal residue
      std::string endICode;    // insertion code of the terminal residue
      int     helixClass;      // helix class
      std::string comment;     // comment about the helix
      int     length;          // length of the helix
      helix_t(int serNum, const std::string &helixID, const std::string &initResName, const std::string &initChainID,
              int initSeqNum, const std::string &initICode,
              const std::string &endResName, const std::string &endChainID, int endSeqNum, const std::string &endICode,
              int helixClass, const std::string &comment, int length) :
         serNum(serNum),
         helixID(helixID),
         initResName(initResName),
         initChainID(initChainID),
         initSeqNum(initSeqNum),
         initICode(initICode),
         endResName(endResName),
         endChainID(endChainID),
         endSeqNum(endSeqNum),
         endICode(endICode),
         helixClass(helixClass),
         comment(comment),
         length(length) {}
   };

   class header_info_t {
   public:
      std::string title;
      std::vector<std::string> journal_lines;
      std::vector<std::string> author_lines;
      std::vector<std::string> compound_lines;
      std::vector<helix_t> helix_info;
      header_info_t() {}
   };
}

#endif // COOT_HEADER_INFO_HH
