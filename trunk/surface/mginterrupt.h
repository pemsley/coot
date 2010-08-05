//   CCP4 Molecular Graphics Program
//
//   Copyright (C) 2004-2005 Stuart McNicholas and  Liz Potterton

//     This code is distributed under the terms and conditions of the
//     CCP4 Program Suite Licence Agreement as a CCP4 Library.
//     A copy of the CCP4 licence can be obtained by writing to the
//     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.

 
 
#ifndef _CCP4MG_INTERRUPT_
#define _CCP4MG_INTERRUPT_

class mginterrupt {
 private:
  static int interrupt_status;
 public:
  mginterrupt(){};
  static void SetStatus(int status_in);
  static int GetStatus(void);
};
#endif
