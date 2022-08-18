//
//  Callbacks.h
//  Aesop
//
//  Created by Martin Noble on 03/04/2012.
//  Copyright (c) 2012 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#ifndef Aesop_Callbacks_h
#define Aesop_Callbacks_h

typedef void (*MyCompletionPtrType)(void *);
typedef void (*MyProgressPtrType)(void *, float);

#ifdef __cplusplus
extern "C" {
#endif
    void globalRedrawCompletionCallback(void *userInfo);   
    void globalRedrawProgressCallback(void *userInfo, float progress);   
    
#ifdef __cplusplus
}
#endif

#endif
