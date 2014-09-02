/*
     mmut/mmut_lipids.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/

#ifndef _CCP4_MMUT_LIPIDS_H_
#define _CCP4_MMUT_LIPIDS_H_
#include <vector>
#include <utility>
#include <mman_manager.h>
#include "atom_util.h"
#include <mgtree.h>

class MMUTLipid{
  /* We have vector of heads here, but currently MMUTLipidCalculate can only cope with 1 head */
  int main_selhnd;
  std::vector<std::vector<Cartesian> > sorted_tail_carts;
  std::vector<std::vector<Cartesian> > head_carts;
  std::vector<std::vector<int> > head_serNums;
  std::vector<std::vector<int> > sorted_tail_serNums;
  std::vector<std::vector<int> > h_tail_serNums;
 public:
  MMUTLipid(){};
  void AddTailCartesians(const std::vector<Cartesian> &carts){sorted_tail_carts.push_back(carts);};
  void AddTailSerNums(const std::vector<int> &carts){sorted_tail_serNums.push_back(carts);};
  void AddTailHSerNums(const std::vector<int> &carts){h_tail_serNums.push_back(carts);};
  void SetHeadCartesians(const std::vector<std::vector<Cartesian> > &head_carts_in){head_carts = head_carts_in;};
  void SetHeadSerNums(const std::vector<std::vector<int> > &head_serNums_in){head_serNums= head_serNums_in;};
  std::vector<std::vector<Cartesian> > GetTailCartesians() const {return sorted_tail_carts;};
  std::vector<std::vector<Cartesian> > GetHeadCartesians() const {return head_carts;};
  std::vector<Cartesian> GetTailCartesian(const int i) const {return sorted_tail_carts[i];};
  std::vector<Cartesian> GetHeadCartesian(const int i) const {return head_carts[i];};
  std::vector<std::vector<int> > GetHeadSerNums() const {return head_serNums;};
  std::vector<int> GetHeadSerNums(const int i) const {return head_serNums[i];};
  std::vector<std::vector<int> > GetTailSerNums() const {return sorted_tail_serNums;};
  std::vector<int> GetTailSerNums(const int i) const {return sorted_tail_serNums[i];};
  std::vector<std::vector<int> > GetTailHSerNums() const {return h_tail_serNums;};
  std::vector<int> GetTailHSerNums(const int i) const {return h_tail_serNums[i];};
  int GetMainSelectionHandle() const {return main_selhnd;};
  void SetMainSelectionHandle(const int i) {main_selhnd=i;};
};

std::vector<MMUTLipid> MMUTLipidCalculate(CMMANManager *molHnd, int selHnd, int minimum_chain_length=8);
int MMUTLipidAnalyse(CMMANManager *molHnd, int selHnd, int minimum_chain_length=8);

#endif //_CCP4_MMUT_LIPIDS_H_
