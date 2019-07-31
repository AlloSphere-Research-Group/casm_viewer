#ifndef AL_VASPREADER
#define AL_VASPREADER

#include <map>
#include <vector>
#include <string>
#include <mutex>

#include "al/core/math/al_Mat.hpp"
#include "al/core/math/al_Vec.hpp"

namespace al {

class VASPReader {
public:
    typedef enum {
        USE_INLINE_ELEMENT_NAMES,   // Set positions using inline names (the name accompanying the position is used instead of following the species declaration line)
        DONT_VALIDATE_INLINE_NAMES       //
    } VASPOption;

    // Made public to allow access to them
    double mNorm {0.0};
    double maxX = 0, maxY = 0, maxZ = 0;
    double minX = 0, minY = 0, minZ = 0;

    VASPReader(std::string basePath = std::string());

    void ignoreElements(std::vector<std::string> elementsToIgnore);

    void setBasePath(std::string path);

    bool loadFile(std::string fileName);

    std::map<std::string, std::vector<float>> &getAllPositions();

    bool hasElement(std::string elementType);

    std::vector<float> &getElementPositions(std::string elementType);

    void setOption(VASPOption option, bool enable = true);

    void print();

    al::Vec3d getNormalizingVector();
    al::Vec3d getCenteringVector();

private:
    std::string mBasePath;
    std::string mFileName;
    bool mVerbose{false};

    std::map<std::string, std::vector<float>> mPositions;
    std::vector<std::string> mElementsToIgnore;

    al::Mat3d mTransformMatrix;
    std::mutex mDataLock;
    std::vector<VASPOption> mOptions;

};

} // namespace al

#endif // AL_VASPREADER
