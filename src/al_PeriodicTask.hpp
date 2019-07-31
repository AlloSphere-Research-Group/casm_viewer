#ifndef AL_PERIODICTASK_HPP
#define AL_PERIODICTASK_HPP

using namespace std;

class PeriodicTask {
public:

  ~PeriodicTask() { stop(); }

  bool start(std::function<bool()> func) {
    if (mFuncThread) {
      return false;
    }
    mRunning = true;
    mFuncThread = std::make_shared<std::thread>(threadFunction,func, this);
    return mFuncThread != nullptr;
  }

  bool running() { return  mRunning || mFuncThread != nullptr; }

  void stop() {
    if (mFuncThread) {
      mRunning = false;
      mFuncThread->join();
      mFuncThread = nullptr;
    }
  }

  std::chrono::nanoseconds waitTime() const;
  void setWaitTime(std::chrono::nanoseconds waitTime);

//  std::chrono::nanoseconds granularity() const
//  {
//    return mGranularity;
//  }
//  void setGranularity(std::chrono::nanoseconds granularity)
//  {
//    mGranularity = granularity;
//  }

private:

  static void threadFunction(std::function<bool()> func, PeriodicTask *thisPtr) {
    auto waitTime = thisPtr->waitTime();
    while(thisPtr->mRunning) {
      if (!func()) {
        thisPtr->mRunning = false;
        break;
      }
      std::this_thread::sleep_for(waitTime);
    }
  }

  std::chrono::nanoseconds mWaitTime {10000ns};
  std::chrono::nanoseconds mGranularity {10000ns}; // Fixme implement granularity for snappier shutdown of waiting
  std::atomic<bool> mRunning;
  std::shared_ptr<std::thread> mFuncThread;

};


std::chrono::nanoseconds PeriodicTask::waitTime() const
{
  return mWaitTime;
}

void PeriodicTask::setWaitTime(std::chrono::nanoseconds waitTime)
{
  mWaitTime = waitTime;
}

#endif // AL_PERIODICTASK_HPP
