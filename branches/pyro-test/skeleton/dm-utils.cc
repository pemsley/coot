
#include <iostream>

using namespace std; 


// fortran indexing.  The (key and ptr) arrays need to be of size
// n+1, since we use key[n].
//
// The sorted postions are 1..n.
// 
void 
shsorti(float *key, int *ptr, int n) {

   int l=1;
   int i;
   int i1, i2;
   int p1, p2; 

   while (l<n) l *=2;
   cout << "l set to " << l << endl;
      
   while (l > 1) {
      l /= 2;
      for (i=1; i<=n-l; i++) {
	 for (i1=i; i1>=1; i1--) {
	    i2 = i1 + 1;
	    p1 = ptr[i1];
	    p2 = ptr[i2];

	    if (key[p2] >= key[p1]) break;
	    
	    cout << "assigning ptr [" << i1 << "] as " << p2 << endl;
	    cout << "assigning ptr [" << i2 << "] as " << p1 << endl;
	    ptr[i1] = p2;
	    ptr[i2] = p1;
	 }
      }
   }
}



