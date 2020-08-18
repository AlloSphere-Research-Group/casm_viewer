#include "modalsynth.hpp"

void ModalVoice::init() {
  residualEnv.levels(0.0f, 1.0f, 0.0f);
  residualEnv.lengths(0.002f, 0.07f);
  setFrequencies(smallHandBell, true);
}

void ModalVoice::onProcess(AudioIOData &io) {
  while (io()) {
    auto ampIt = amps.begin();
    float excitation = noise() * residualEnv();
    for (auto &mode : modes) {
      float out = mode(excitation);
      io.out(0) += *ampIt++ * globalAmp * out;
    }
    envFollow(io.out(0));
  }
  if (envFollow.done(0.0001)) {
    free();
  }
}

void ModalVoice::onTriggerOn() {
  residualEnv.reset();

  for (auto &mode : modes) {
    mode.zero();
  }
}

void ModalVoice::setFrequencies(std::vector<float> freqs,
                                float autoWidthFactor) {
  modes.resize(freqs.size());
  amps.resize(freqs.size());
  auto modesIt = modes.begin();
  auto ampsIt = amps.begin();
  int counter = 1;
  for (auto f : freqs) {
    modesIt->freq(f * fundamentalFreq);
    //      xylo 0.006
    // aluminium 0.00012
    // tubularBell 0.0004
    // small handbell 0.0004
    // small handbell 0.005 // low frequencies sounds like a pot
    if (autoWidthFactor > 0) {
      modesIt->width(f * fundamentalFreq * autoWidthFactor);
    }
    *ampsIt = 1.0 / counter++;
    modesIt++;
    ampsIt++;
  }
}
