#ifndef PARAMETERSPACE_HPP
#define PARAMETERSPACE_HPP

#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#endif

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "al/ui/al_Parameter.hpp"

namespace al {

/**
 * @brief The ParameterSpace class maps parameter values to string ids
 *
 * This allows mapping continuous parameters to string ids, for example
 * for mapping to directory structures
 *
 * ParameterSpaces can be linked together using addConnectedParameterSpace()
 * to have the final ids be deduces from the combinations of two parameters
 *
 */
class ParameterSpace {
public:
  ParameterSpace(std::string name, std::string group = "")
      : mParameterValue(name, group) {}

  // Access to current

  float getCurrentValue();

  void setCurrentValue(float value);

  size_t getCurrentIndex();

  std::string getCurrentId();

  // the parameter instance holds the current value.
  // You can set values for parameter space through this function
  // Register notifications and create GUIs/ network synchronization
  // Through this instance.
  Parameter &parameter();

  // Multidimensional parameter spaces will result in single values having
  // multiple ids. This can be resolved externally using this function
  // but perhaps we should have a higher level class that solves this issue for
  // the general case
  std::vector<std::string> getAllCurrentIds();

  std::vector<size_t> getAllCurrentIndeces();

  // Access to specific elements

  // When using reverse = true, the value returned describes the
  // end of an open range, i.e, the index returned is one more than
  // an index that would match the value.
  size_t getFirstIndexForValue(float value, bool reverse = false);

  float at(size_t x);

  std::string idAt(size_t x);

  size_t getIndexForValue(float value);

  // Access to complete sets

  std::vector<std::string> getAllIds(float value);

  std::vector<size_t> getAllIndeces(float value);

  std::vector<std::pair<std::string, float>> values();

  // Move parameter space

  void stepIncrement();

  void stepDecrease();

  //
  size_t size();

  // Modification of parameter space. These are not protected by the lock, the
  // user needs to lock as needed
  void sort(std::function<bool(const std::pair<std::string, float> &a,
                               const std::pair<std::string, float> &b)>
                sortFunction =
                    [](const std::pair<std::string, float> &a,
                       const std::pair<std::string, float> &b) -> bool {
    return a.second < b.second;
  });

  void clear();

  void push_back(float value, std::string id = "");

  void addConnectedParameterSpace(ParameterSpace *paramSpace);

  // Protect parameter space (to avoid access during modification)

  // TODO all readers above need to use this lock
  void lock() { mLock.lock(); }

  void unlock() { mLock.unlock(); }

private:
  // Data
  std::vector<std::pair<std::string, float>> mValues;

  std::mutex mLock;

  // Current state
  Parameter mParameterValue;

  std::vector<ParameterSpace *> mConnectedSpaces;
};

} // namespace al

#endif // PARAMETERSPACE_HPP
