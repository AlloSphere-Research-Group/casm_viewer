#include "parameterspace.hpp"

using namespace al;

float ParameterSpace::at(size_t x) {
  if(x < mValues.size()) {
    auto it = mValues.begin();
    std::advance(it,x);
    return it->second;
  } else {
    return 0.0f;
  }
}

std::string ParameterSpace::idAt(size_t x) {
  if(x < mValues.size()) {
    auto it = mValues.begin();
    std::advance(it,x);
    return it->first;
  } else {
    return std::string();
  }
}

size_t ParameterSpace::size() {
  return mValues.size();
}

void ParameterSpace::clear() {
  mParameterValue.min(FLT_MAX);
  mParameterValue.max(FLT_MIN);
  mValues.clear();
}

size_t ParameterSpace::getIndexForValue(float value) {
  auto indeces = getAllIndeces(value);
  for (auto connectedSpace: mConnectedSpaces) {
    if (connectedSpace->size() > 0) {
      std::vector<std::string> ids = connectedSpace->getAllIds(connectedSpace->parameter().get());
      std::vector<size_t> commonIndeces;
      for (auto index: indeces) {
        auto id = idAt(index);
        if (std::find(ids.begin(), ids.end(), idAt(index)) != ids.end()) {
          commonIndeces.push_back(index);
        } else if (ids.size() == 1 && ids[0] == "./") {
          commonIndeces.push_back(index);
        }
      }
      indeces = commonIndeces;
    }
  }
  if (indeces.size() == 1) {
    return indeces[0];
  }
  return -1;
}

size_t ParameterSpace::getFirstIndexForValue(float value, bool reverse) {
  int paramIndex = -1;

  if (!reverse) {
    size_t i = 0;
    for (auto it = mValues.begin(); it != mValues.end(); it++) {
      if (it->second == value) {
        paramIndex = i;
        break;
      } else if (it->second > value
                 && (i == mValues.size() - 1)) {
        break;
      }
      auto next = it;
      std::advance(next, 1);
      if (next != mValues.end()) {
        if (it->second > value && next->second < value) { // space is sorted and descending
          paramIndex = i;
          break;
        } else if (it->second < value && next->second > value) { // space is sorted and ascending
          paramIndex = i + 1;
          break;
        }

      }
      i++;
    }
  } else {
    value = at(getFirstIndexForValue(value));
    size_t i = mValues.size();
    for (auto it = mValues.rbegin(); it != mValues.rend(); it++) {
      if (it->second == value) {
        paramIndex = i;
        break;
      } else if (it->second < value
                 && (i == 0)) {
        break;
      }
      auto next = it;
      std::advance(next, 1);
      if (next != mValues.rend()) {
        if (it->second < value && next->second > value) { // space is sorted and descending
          paramIndex = i;
          break;
        } else if (it->second > value && next->second < value) { // space is sorted and ascending
          paramIndex = i - 1;
          break;
        }

      }
      i--;
    }
  }
  if (paramIndex < 0) {
    //          std::cerr << "WARNING: index not found" << std::endl;
    paramIndex = 0;
  }
  return paramIndex;
}

float ParameterSpace::getCurrentValue() {
  return mValues[getCurrentIndex()].second;
}

void ParameterSpace::setCurrentValue(float value) {
  parameter().set(value);
}

size_t ParameterSpace::getCurrentIndex() {
  return getIndexForValue(mParameterValue.get());
}

std::string ParameterSpace::getCurrentId() {
  return idAt(getCurrentIndex());
}

std::vector<std::string> ParameterSpace::getAllCurrentIds() {
  float value = getCurrentValue();
  return getAllIds(value);
}

std::vector<std::string> ParameterSpace::getAllIds(float value) {
  size_t lowIndex = getFirstIndexForValue(value);
  size_t highIndex = getFirstIndexForValue(value, true); // Open range value (excluded from range)
  std::vector<std::string> ids;

  for(size_t i = lowIndex; i < highIndex; i++) {
    ids.push_back(idAt(i));
  }
  if (lowIndex == highIndex) { // Hack... shouldn't need to be done. This should be fixed for one dimensional data parameter spaces.
    ids.push_back(idAt(lowIndex));
  }
  return ids;
}

std::vector<size_t> ParameterSpace::getAllCurrentIndeces() {
  float value = getCurrentValue();
  return getAllIndeces(value);
}

std::vector<size_t> ParameterSpace::getAllIndeces(float value) {
  size_t lowIndex = getFirstIndexForValue(value);
  size_t highIndex = getFirstIndexForValue(value, true); // Open range value (excluded from range)
  std::vector<size_t> idxs;

  for(size_t i = lowIndex; i < highIndex; i++) {
    idxs.push_back(i);
  }
  if (lowIndex == highIndex) { // Hack... shouldn't need to be done. This should be fixed for one dimensional data parameter spaces.
    idxs.push_back(lowIndex);
  }
  return idxs;
}

Parameter &ParameterSpace::parameter() {
  return mParameterValue;
}

void ParameterSpace::stepIncrement() {
  size_t curIndex = getCurrentIndex();
  float temp = mParameterValue.get();
  if (int64_t(curIndex) < int64_t(mValues.size()) - 1) { // Check if we have an element above to compare to
    size_t offset = 1;
    float nextTemp = at(curIndex + offset);
    while (nextTemp == temp && curIndex < mValues.size() - ++offset) {
      nextTemp = at(curIndex + offset);
    }
    if (nextTemp > temp) {// Element above is greater, so increment index and load
      parameter().set(nextTemp);
    } else if (curIndex > 0 && nextTemp != temp) {// If element above is not greater and we have space to go down, then decrement index and load
      float prevTemp = at(curIndex - 1);
      parameter().set(prevTemp);
    }
  } else if (curIndex > 0) { // Can we move one below?
    float prevTemp = at(curIndex - 1);
    if (prevTemp > temp) {
      parameter().set(prevTemp);
    }
  }
}

void ParameterSpace::stepDecrease() {
  size_t curIndex = getCurrentIndex();
  float temp = mParameterValue.get();
  float nextTemp = temp;
  while (nextTemp == temp  && (int) curIndex < (int) mValues.size() - 1) {
    nextTemp = at(++curIndex);
  }
  if (nextTemp < temp) {// Element above is greater, so increment index and load
    parameter().set(nextTemp);
  } else {
    float prevTemp = temp;
    while (prevTemp == temp  && (int) curIndex > 1) {
      prevTemp = at(--curIndex);
    }
    if (prevTemp != temp) {
      parameter().set(prevTemp);
    }
  }
}

void ParameterSpace::push_back(float value, std::string id) {
  if (id == "") {
    std::stringstream ss; // Use stringstream for better handling of trailing zeros
    ss << value;
    id = ss.str();
  }

  bool found = false;
  for (auto it = mValues.begin(); it != mValues.end(); it++) {
    if (it->first == id) {
      it->second = value;
      found = true;
      break; // Assume ids are unique
    }
  }

  if (!found) {
    mValues.emplace_back(std::pair<std::string, float>(id, value));
  }

  if (value > mParameterValue.max()) {
    mParameterValue.max(value);
  }
  if (value < mParameterValue.min()) {
    mParameterValue.min(value);
  }
}

void ParameterSpace::addConnectedParameterSpace(ParameterSpace *paramSpace) {
  mConnectedSpaces.push_back(paramSpace);
}

void ParameterSpace::sort(std::function<bool(const std::pair<std::string, float> & a, const std::pair<std::string, float> & b)> sortFunction) {
  std::sort(mValues.begin(), mValues.end(),  sortFunction);
}

std::vector<std::pair<std::string, float> > ParameterSpace::values() {
  return mValues;
}
