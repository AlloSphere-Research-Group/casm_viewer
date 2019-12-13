
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>

#include "al/io/al_File.hpp"

#include "al_VASPReader.hpp"

#ifdef AL_WINDOWS
#include <Windows.h>
#endif

using namespace al;

bool VASPReader::loadFile(std::string fileName) {
  std::unique_lock<std::mutex> lk(mDataLock);
  std::string fullName = al::File::conformPathToOS(mBasePath);
  if (mBasePath == "./") {
    fullName = "";
  }
  fullName += fileName;
  std::ifstream infile(fullName);
  if (infile.fail()) {
    std::cout << "ERROR: Can't access file " << fullName << std::endl;
    return false;
  }
  mFileName = fileName;
  //        std::cout << "Reading file: " << fullName << std::endl;

  std::istringstream ss;

  std::string line;
  std::getline(infile, line);
  ss.str(line);
  std::string format;
  if (line.size() > 0 && !(ss >> format)) {
    ss.clear();
    std::cout << "Error. expecting FORMAT instead of " << line << std::endl;
  }

  std::getline(infile, line);
  ss.str(line);
  double scale;
  if (!(ss >> scale)) {
    ss.clear();
    std::cout << "Error. expecting SCALE instead of " << line << std::endl;
  }

  // Matrix
  std::getline(infile, line);
  ss.str(line);
  ss.clear();

  if (!(ss >> mTransformMatrix(0, 0) >> mTransformMatrix(0, 1) >>
        mTransformMatrix(0, 2))) {
    ss.clear();
    std::cerr << "Error. expecting MATRIXROW";
    if (mVerbose) {
      std::cout << "instead of " << line;
    }
    std::cerr << std::endl;
  }

  std::getline(infile, line);
  ss.str(line);
  ss.clear();
  if (!(ss >> mTransformMatrix(1, 0) >> mTransformMatrix(1, 1) >>
        mTransformMatrix(1, 2))) {
    ss.clear();
    std::cerr << "Error. expecting MATRIXROW";
    if (mVerbose) {
      std::cout << "instead of " << line;
    }
    std::cerr << std::endl;
  }

  std::getline(infile, line);
  ss.str(line);
  ss.clear();
  if (!(ss >> mTransformMatrix(2, 0) >> mTransformMatrix(2, 1) >>
        mTransformMatrix(2, 2))) {
    ss.clear();
    std::cerr << "Error. expecting MATRIXROW";
    if (mVerbose) {
      std::cout << "instead of " << line;
    }
    std::cerr << std::endl;
  }

  // Element Names
  std::getline(infile, line);
  ss.str(line);
  ss.clear();
  std::string buf;
  std::vector<std::string> elementNames;
  while (ss >> buf) {
    elementNames.push_back(buf);
  }
  if (elementNames.size() == 0) {
    std::cout << "ERROR: No elements listed in VASP file" << std::endl;
    return false;
  }

  ss.clear();
  //        mTransformMatrix *= 0.001/(49*2.0);

  // Element Counts
  std::getline(infile, line);
  ss.str(line);
  int bufint;
  std::vector<int> elementCount;
  //        int i = 0;
  while (ss >> bufint) {
    elementCount.push_back(bufint);
  }

  //        for (auto e: elementCount) {
  //            std::cout << e.first << "..." << e.second << std::endl;
  //        }

  std::getline(infile, line);
  std::string mode = line;
  bool useTransformMatrix;
  if (mode == "Direct") {
    useTransformMatrix = true;
  } else if (mode == "Cartesian") {
    useTransformMatrix = false;
  } else {
    std::cout << "VASPReader: Unrecognized mode " << mode
              << ". Assuming Cartesian." << std::endl;
    useTransformMatrix = false;
  }

  std::map<std::string, unsigned int> counter;

  mPositions.clear();
  for (unsigned int i = 0; i < elementNames.size(); i++) {
    auto &name = elementNames[i];
    auto &count = elementCount[i];
    if (std::find(mElementsToIgnore.begin(), mElementsToIgnore.end(), name) ==
        mElementsToIgnore.end()) {
      if (mPositions.find(name) == mPositions.end()) {
        mPositions[name] = std::vector<float>();
      }
      mPositions[name].resize(mPositions[name].size() +
                              count * 4);  // Pre-allocate for faster reading
      counter[name] = 0;
    }
  }

  mNorm = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      mNorm += mTransformMatrix(i, j) * mTransformMatrix(i, j);
    }
  }
  mNorm = sqrt(mNorm);

  maxX = 0.0;
  maxY = 0.0;
  maxZ = 0.0;
  minX = FLT_MAX;
  minY = FLT_MAX;
  minZ = FLT_MAX;

  bool useInlineNames = std::find(mOptions.begin(), mOptions.end(),
                                  USE_INLINE_ELEMENT_NAMES) != mOptions.end();
  bool validateInlineNames =
      std::find(mOptions.begin(), mOptions.end(), DONT_VALIDATE_INLINE_NAMES) ==
      mOptions.end();
  auto currentTypeIt = elementNames.begin();
  auto currentTypeNumberIt = elementCount.begin();
  int currentElementCount = *currentTypeNumberIt;
  std::string currentType = *currentTypeIt;
  int totalCounter = 0;
  std::cout << "Parsing atoms" << std::endl;
  while (std::getline(infile, line)) {
    al::Vec3d pos;

    std::string type;
    ss.str(line);
    ss.clear();
    if (line.size() > 0) {
      if (!(ss >> pos.x >> pos.y >> pos.z >> type)) {
        ss.clear();
        std::cout << "Error, unexpected: " << line << std::endl;
      } else {
        //                    pos /= mNorm;
        if (useInlineNames) {
          currentType = type;
          if (mPositions[currentType].size() >= counter[currentType] * 4) {
            mPositions[currentType].reserve(mPositions[currentType].size() +
                                            1024);
            mPositions[currentType].resize(mPositions[currentType].size() + 4);
          } else {
            std::cout << "POSCAR size mismatch." << std::endl;
          }
        } else {
          if (totalCounter >= currentElementCount) {
            currentTypeIt++;
            currentTypeNumberIt++;
            if (currentTypeIt != elementNames.end() ||
                currentTypeNumberIt != elementCount.end()) {
              currentType = *currentTypeIt;
              currentElementCount += *currentTypeNumberIt;
            } else {
              currentTypeIt--;
              currentTypeNumberIt--;
              std::cerr << "POSCAR size mismatch" << std::endl;
            }
            //                            counter[currentType] = 0;
          }
        }
        if (validateInlineNames) {
          if (currentType != type && type != "?" && type != "X") {
            std::cout << "VASPReader: species mismatch. Got " << type
                      << " expecting " << currentType << std::endl;
          }
        }
        if (std::find(mElementsToIgnore.begin(), mElementsToIgnore.end(),
                      currentType) == mElementsToIgnore.end()) {
          double x;
          if (useTransformMatrix) {
            x = pos.dot(mTransformMatrix.col(0));
          } else {
            x = pos.x;
          }
          assert(mPositions[currentType].size() > counter[currentType] * 4 + 3);
          mPositions[currentType][counter[currentType] * 4] = x;
          if (maxX < x) {
            maxX = x;
          }
          if (minX > mPositions[currentType][counter[currentType] * 4]) {
            minX = mPositions[currentType][counter[currentType] * 4];
          }
          double y;
          if (useTransformMatrix) {
            y = pos.dot(mTransformMatrix.col(1));
          } else {
            y = pos.y;
          }
          mPositions[currentType][counter[currentType] * 4 + 1] = y;
          if (maxY < y) {
            maxY = y;
          }
          if (minY > y) {
            minY = y;
          }
          double z;
          if (useTransformMatrix) {
            z = pos.dot(mTransformMatrix.col(2));
          } else {
            z = pos.z;
          }
          mPositions[currentType][counter[currentType] * 4 + 2] = z;
          if (maxZ < z) {
            maxZ = z;
          }
          if (minZ > z) {
            minZ = z;
          }
          mPositions[currentType][counter[currentType] * 4 + 3] = 0.0;
        }
        counter[currentType]++;
        totalCounter++;
      }
    }
  }

  std::cout << "Done Parsing atoms" << std::endl;
  //        for (auto counted: counter) {
  //            std::cout << "VASPReader Read " << counted.first << ":" <<
  //            counted.second <<std::endl;
  //        }
  infile.close();
  return !infile.bad();
}

std::vector<float> &VASPReader::getElementPositions(std::string elementType) {
  std::unique_lock<std::mutex> lk(mDataLock);
  if (mPositions.find(elementType) == mPositions.end()) {
    mPositions[elementType] = std::vector<float>();
    std::cout << "ERROR: Invalid element: " << elementType << std::endl;
  }
  return mPositions[elementType];
}

void VASPReader::setOption(VASPReader::VASPOption option, bool enable) {
  if (enable) {
    if (find(mOptions.begin(), mOptions.end(), option) == mOptions.end()) {
      mOptions.push_back(option);
    }
  } else {  // remove option from option list
    if (find(mOptions.begin(), mOptions.end(), option) == mOptions.end()) {
      std::remove(mOptions.begin(), mOptions.end(), option);
    }
  }
}

void VASPReader::print() {
  std::cout << "Path: " << mBasePath << std::endl;
  std::cout << "File: " << mFileName << std::endl;
  std::cout << "Elements: ";
  for (auto position : mPositions) {
    std::cout << position.first << ":" << position.second.size() / 4 << " ";
  }
  std::cout << "Matrix: ";

  std::cout << "Max: " << maxX << "," << maxY << "," << maxZ << std::endl;
  std::cout << "Min: " << minX << "," << minY << "," << minZ << std::endl;
  std::cout << std::endl;
}

VASPReader::VASPReader(std::string basePath) { setBasePath(basePath); }

void VASPReader::ignoreElements(std::vector<std::string> elementsToIgnore) {
  mElementsToIgnore = elementsToIgnore;
}

void VASPReader::setBasePath(std::string path) { mBasePath = path; }

std::map<std::string, std::vector<float> > &VASPReader::getAllPositions() {
  return mPositions;
}

bool VASPReader::hasElement(std::string elementType) {
  return mPositions.find(elementType) != mPositions.end();
}

al::Vec3d VASPReader::getNormalizingVector() {
  return al::Vec3d(0.9 / (maxX - minX), 0.9 / (maxY - minY),
                   0.9 / (maxZ - minZ));
}

al::Vec3d VASPReader::getCenteringVector() {
  return al::Vec3d((maxX + minX) / 2.0, (maxY + minY) / 2.0,
                   (maxZ + minZ) / 2.0) /
         mNorm;
}
