#include <vector>
// helper trait to detect std::vector
template <typename T>
struct is_vector : std::false_type {};

template <typename T, typename A>
struct is_vector<std::vector<T, A>> : std::true_type {};


// main recursive mapf
template <typename Container, typename F>
auto mapf(const Container& source, F&& function) {
  if constexpr (is_vector<Container>::value) {
    using InnerType = typename Container::value_type;
    std::vector<decltype(mapf(std::declval<InnerType>(), std::forward<F>(function)))> result;
    result.reserve(source.size());
    for (auto&& elem : source)
      result.push_back(mapf(elem, std::forward<F>(function)));
    return result;
  } else {
    // base case: not a vector, apply the function directly
    return function(source);
  }
}
