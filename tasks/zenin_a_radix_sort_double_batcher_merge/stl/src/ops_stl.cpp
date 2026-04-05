#include "zenin_a_radix_sort_double_batcher_merge/stl/include/ops_stl.hpp"

//#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

//#include "oneapi/tbb/parallel_for.h"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

ZeninARadixSortDoubleBatcherMergeSTL::ZeninARadixSortDoubleBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeSTL::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSTL::PreProcessingImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSTL::RunImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge