
namespace Filter {

template <typename NumType>
struct FilterStructure {
  unsigned int kernel_size;
  unsigned int os_buf_len;
  unsigned int circular_buf_len;
  double stopband_amp_reduction;
  float transition_bandwidth;
  NumType *os_buf;
  NumType *kernel;
  NumType *circular_buf;
};

} // namespace Filter 

