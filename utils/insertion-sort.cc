
#include <algorithm>

#include "insertion-sort.hh"

void insertionSort(std::vector<int> &vec) {
   for (auto it = vec.begin(); it != vec.end(); it++) {
      // Searching the upper bound, i.e., first element greater than *it from beginning
      auto const insertion_point = std::upper_bound(vec.begin(), it, *it);
      std::rotate(insertion_point, it, it+1);
   }
}
