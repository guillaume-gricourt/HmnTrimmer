// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Define trimmers
// ============================================================================
#ifndef APP_HMNTRIMMER_TRIMMER_H_
#define APP_HMNTRIMMER_TRIMMER_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <cstring>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "seqs.hpp"

using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag QualTail
// ----------------------------------------------------------------------------

struct QualTail_;
typedef Tag<QualTail_> QualTail;

// ----------------------------------------------------------------------------
// Tag QualSld
// ----------------------------------------------------------------------------

struct QualSld_;
typedef Tag<QualSld_> QualSld;

// ----------------------------------------------------------------------------
// Tag LenMin
// ----------------------------------------------------------------------------

struct LenMin_;
typedef Tag<LenMin_> LenMin;

// ----------------------------------------------------------------------------
// Tag InfoDust
// ----------------------------------------------------------------------------

struct InfoDust_;
typedef Tag<InfoDust_> InfoDust;

// ----------------------------------------------------------------------------
// Tag InfoN
// ----------------------------------------------------------------------------

struct InfoN_;
typedef Tag<InfoN_> InfoN;

// ----------------------------------------------------------------------------
// TagList Trimmers
// ----------------------------------------------------------------------------

typedef TagList<QualTail,
    TagList<QualSld,
    TagList<LenMin,
    TagList<InfoDust,
    TagList<InfoN
    > > > > > Trimmers;

// ============================================================================
// Classes
// ============================================================================

template <typename TTrimmer, typename T = void>
struct IdTrimmer;

// ----------------------------------------------------------------------------
// Class QualTail
// ----------------------------------------------------------------------------

template <typename T>
struct IdTrimmer<QualTail_, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * IdTrimmer<QualTail_, T>::VALUE[1] =
{
    "QualTail"
};

// ----------------------------------------------------------------------------
// Class QualSld
// ----------------------------------------------------------------------------

template <typename T>
struct IdTrimmer<QualSld_, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * IdTrimmer<QualSld_, T>::VALUE[1] =
{
    "QualSld"
};

// ----------------------------------------------------------------------------
// Class LenMin
// ----------------------------------------------------------------------------

template <typename T>
struct IdTrimmer<LenMin_, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * IdTrimmer<LenMin_, T>::VALUE[1] =
{
    "LenMin"
};

// ----------------------------------------------------------------------------
// Class InfoDust
// ----------------------------------------------------------------------------

template <typename T>
struct IdTrimmer<InfoDust_, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * IdTrimmer<InfoDust_, T>::VALUE[1] =
{
    "InfoDust"
};

// ----------------------------------------------------------------------------
// Class InfoN
// ----------------------------------------------------------------------------

template <typename T>
struct IdTrimmer<InfoN_, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * IdTrimmer<InfoN_, T>::VALUE[1] =
{
    "InfoN"
};

// ----------------------------------------------------------------------------
// Class InfoValues
// ----------------------------------------------------------------------------

template <typename TValue = unsigned, typename TFValue = float>
struct InfoValues
{

    static inline TValue const getWindowSize()
    {
        TValue const _data = 64;
        return _data;
    }
    static inline TValue const getWindowStep()
    {
        TValue const _data = 32;
        return _data;
    }

    static inline TFValue const getWindowMax()
    {
        TFValue const _data = 62.0;
        return _data;
    }

    static inline TFValue const getByNum()
    {
        TFValue const _data = 1.0 / getWindowMax();
        return _data;
    }

    static inline TFValue const getOneOverLog()
    {
        TFValue const _data = 1.0 / std::log(62);
        return  _data;
    }
};

#endif  // #ifndef APP_HMNTRIMMER_TRIMMER_H_
