// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Timer
// ============================================================================
#ifndef APP_HMNTRIMMER_TIMER_H_
#define APP_HMNTRIMMER_TIMER_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <chrono>
#include <iostream>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Timer
// ----------------------------------------------------------------------------

template <typename TClock = std::chrono::high_resolution_clock>
struct Timer
{
    typedef std::chrono::time_point<TClock> TTime;

    TTime begin;
    TTime end;

    Timer() :
        begin(),
        end()
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function start()
// ----------------------------------------------------------------------------

template <typename TClock>
inline void start(Timer<TClock> & timer)
{
    timer.begin = TClock::now();
}

// ----------------------------------------------------------------------------
// Function stop()
// ----------------------------------------------------------------------------

template <typename TClock>
inline void stop(Timer<TClock> & timer)
{
    timer.end = TClock::now();
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TClock, typename TValue>
inline void getTimer(Timer<TClock> const & timer, TValue & time)
{
    time = std::chrono::duration_cast<std::chrono::seconds>(timer.end - timer.begin).count();
}

#endif // APP_HMNTRIMMER_TIMER_H_
