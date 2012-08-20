#ifndef MRR_PDB_UTILITY_HXX_
#define MRR_PDB_UTILITY_HXX_

#include <tuple>
#include <chrono>
#include <string>
#include <type_traits>
#include <regex>

namespace mrr {
namespace util {

template <typename ClockType = std::chrono::high_resolution_clock>
struct stopwatch
{
private:
  using time_point = typename ClockType::time_point;

public:
  stopwatch()
    : total_time_span(0)
  {
  }

  double lap()
  {
    using std::chrono::duration;
    using std::chrono::duration_cast;

    time_point lap_end_time = ClockType::now();
    lap_time_span = duration_cast<duration<double> >(
      lap_end_time - lap_start_time).count();
    lap_start_time = lap_end_time;
    return lap_time_span;
  }

  void reset()
  {
    total_time_span = 0;
  }

  void start()
  {
    start_time = ClockType::now();
    lap_start_time = start_time;
  }

  double stop()
  {
    using std::chrono::duration;
    using std::chrono::duration_cast;

    time_point end_time = ClockType::now();
    total_time_span += duration_cast<duration<double> >(end_time - start_time).count();
    return total_time_span;
  }

private:
  time_point start_time;
  double total_time_span;
  time_point lap_start_time;
  double lap_time_span;
};


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename Clock, typename Func>
inline auto time_function(Func&& f)
  -> typename std::enable_if<
       !std::is_void<decltype(f())>::value,
       std::tuple<
         decltype(f()),
         typename Clock::duration
       >
     >::type
{
    auto start_time = Clock::now();
    auto result = f();
    auto time_taken = Clock::now() - start_time;
    return std::make_tuple(result, time_taken);
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename Clock, typename Func>
inline auto time_function(Func&& f)
  -> typename std::enable_if<
       std::is_void<decltype(f())>::value,
       typename Clock::duration
     >::type
{
    auto start_time = Clock::now();
    f();
    auto time_taken =  Clock::now() - start_time;
    return time_taken;
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename T>
std::istream& ignore(std::istream& is)
{
  T t;
  is >> t;
  return is;
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool has_prefix(std::string const& input, std::string const& prefix)
{
  return (input.find(prefix) == 0);
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
}
}


#endif // #ifndef MRR_PDB_UTILITY_HXX_
