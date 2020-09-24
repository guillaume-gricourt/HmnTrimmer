// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Store sequences and trims
// ============================================================================
#ifndef APP_HMNTRIMMER_SEQS_H_
#define APP_HMNTRIMMER_SEQS_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <set>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/alignment_free.h> 
#include <seqan/misc/accumulators.h> 
#include <seqan/seq_io.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "misc_tags.hpp"
#include "trimmers.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeqConfig
// ----------------------------------------------------------------------------

template <typename TRunConfig>
struct SeqConfig
{
    typedef CharString          TAlphabet;
    typedef Dna5Q               TAlphabetSequence;
    typedef Alloc<>             TSeqSpec;

    typedef Owner<>             TSeqsNameSpec;
    typedef Owner<>             TSeqsSpec;

    typedef typename TRunConfig::TThreading         TThreading;
    typedef typename TRunConfig::TSequencing        TSequencing;
    typedef typename TRunConfig::TInputFormat       TInputFormat;
    typedef typename TRunConfig::TOutputFormat      TOutputFormat;
};

// ----------------------------------------------------------------------------
// Class SeqStore
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct SeqStore
{
    //static constexpr const int ReadRecords = 100000; Using : SeqStore<TSpec,
    //TConfig>::ReadRecords
    typedef typename TConfig::TAlphabet                     TAlphabet;
    typedef typename TConfig::TAlphabetSequence             TAlphabetSequence;    
    typedef typename TConfig::TSeqSpec                      TSeqSpec;

    typedef typename TConfig::TSeqsNameSpec                 TSeqsNameSpec;
    typedef typename TConfig::TSeqsSpec                     TSeqsSpec;

    typedef typename TConfig::TThreading                    TThreading;
    typedef typename TConfig::TSequencing                   TSequencing;
    typedef typename TConfig::TInputFormat                  TInputFormat;
    typedef typename TConfig::TOutputFormat                 TOutputFormat;

    // Global parameters.
    typedef StringSet<TAlphabet, TSeqsNameSpec>             TSeqName;
    typedef String<TAlphabetSequence, TSeqSpec>             TSeq;
    typedef StringSet<TSeq, TSeqsSpec>                      TSeqs;

    // Idents.
    typedef typename Iterator<TSeqs const, Standard>::Type  TSeqsIt;
    typedef typename Size<TSeqs>::Type                      TReadId;
    typedef typename Value<TSeqs const>::Type               TSeqsValue;
    typedef typename Size<TSeqsValue>::Type                 TSize;
    typedef typename Id<TSeqs>::Type                        TId;

    typedef std::set<TId>                                   TIdents;

    typedef Triple<TSeqName>                                TPName;
    typedef Triple<TSeqs>                                   TPSeqs;
    
    TPName      names;
    TPSeqs      seqs;

    TPName      namesRemove;
    TPSeqs      seqsRemove;

    TIdents     idents;

    SeqStore() :
        names(),
        seqs(),
        namesRemove(),
        seqsRemove(),
        idents()
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function trimPrefix()
// ----------------------------------------------------------------------------

template <typename TString, typename TValue>
inline typename Prefix<TString const>::Type
trimPrefix(TString const & string, TValue const & pos)
{
    return prefix(string, pos);
}

// ----------------------------------------------------------------------------
// Function update()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIdents>
inline void update(TContainer & c, TIdents & idents)
{
    // Adapted from https://codereview.stackexchange.com - 206686
    auto it = std::begin(idents);
    unsigned current_index = 0;

    assert(std::is_sorted(it, std::end(idents)));
    typedef typename Value<TContainer const>::Type    TContainerValue;
    typename Iterator<TContainer, Standard>::Type newEnd = std::remove_if(
            begin(c, Standard()), 
            end(c, Standard()), 
            [&](TContainerValue const &){
            if (it == std::end(idents)) { return false; }
            if (*it == current_index++) { return ++it, true; }
            else {return false;}}
    );
    resize(c, position(newEnd, c) , Exact());
}

template <typename TSpec, typename TConfig>
inline void update(SeqStore<TSpec, TConfig> & me)
{
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {    
            update(me.names.i1, me.idents);
        }
        SEQAN_OMP_PRAGMA(section)
        {
            update(me.seqs.i1, me.idents);
        }
    }

    if(IsSameType<typename TConfig::TSequencing, SequencingPaired>::VALUE)
    {
        SEQAN_OMP_PRAGMA(parallel sections)
        {
            SEQAN_OMP_PRAGMA(section)
            {    
                update(me.names.i2, me.idents);
            }
            SEQAN_OMP_PRAGMA(section)
            {
                update(me.seqs.i2, me.idents);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function updateDiscard()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIdents>
inline void
updateDiscard(TContainer & origin, TContainer & incoming, TIdents & idents)
{
    typedef typename TIdents::const_iterator TIdentIt;
    for (TIdentIt itIdent = idents.begin(); itIdent != idents.end(); ++ itIdent)
    {
            appendValue(incoming, origin[*itIdent]);
    }
}


template <typename TSpec, typename TConfig>
inline void updateDiscard(SeqStore<TSpec, TConfig> & me)
{
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {    
            updateDiscard(me.names.i1, me.namesRemove.i1, me.idents);
        }
        SEQAN_OMP_PRAGMA(section)
        {
            updateDiscard(me.seqs.i1, me.seqsRemove.i1, me.idents);
        }
    }

    if(IsSameType<typename TConfig::TSequencing, SequencingPaired>::VALUE)
    {   
        SEQAN_OMP_PRAGMA(parallel sections)
        {
            SEQAN_OMP_PRAGMA(section)
            {    
                updateDiscard(me.names.i2, me.namesRemove.i2, me.idents);
            }
            SEQAN_OMP_PRAGMA(section)
            {
                updateDiscard(me.seqs.i2, me.seqsRemove.i2, me.idents);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function trim() - QualTail
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TIdents, typename TMap, typename TThreading>
inline void 
trim(TSeqs & seqs, TIdents & idents, TMap & map, TThreading const & threading,
QualTail)
{
    typedef typename Iterator<TSeqs const, Standard>::Type      TSeqsIt;
    typedef typename Value<TSeqs const>::Type                   TSeqsValue;
    typedef typename Position<TSeqs>::Type                      TReadPos;
    typedef typename Size<TSeqsValue>::Type                     TSize;
    typedef typename Id<TSeqs>::Type                            TId;
    typedef typename TMap::mapped_type                          TParam;
    
    int baseQuality = static_cast<int>(map.at("base_quality"));
    TParam baseNumber = map.at("base_number");
    int lenPerc = -1;
    if (map.count("len_perc"))
        lenPerc = static_cast<int>(map.at("len_perc"));
    
    // Iterate.
    iterate(seqs, [&](TSeqsIt const & it)
    {
        TReadPos readPos = position(it, seqs);
        TId id = positionToId(seqs, readPos);
        TSeqsValue const & matches = value(it);
        TSize size = length(matches);

        // Discard record.
        if(size < baseNumber)
        {
            idents.insert(id);
            return;
        }

        // Main loop.
        TSize lentokeep = size;
        TSize count = 0;
        for(int i = size-1; i > -1 ; --i)
        {
            if (getQualityValue(matches[i]) <= baseQuality)
            {
                    ++count;
            }
            else
            {
                    count = 0;
            }

            if (count >= baseNumber)
            {
                lentokeep = i;
            }
        }                    

        // Choose keeping or discarding record.
        if (lentokeep != size)
        {
            if(lentokeep < 1 or 
                (lenPerc > 0 and lentokeep*100.0/size < lenPerc))
            {
                SEQAN_OMP_PRAGMA(critical)
                {
                    idents.insert(id);
                }
            }
            else
            {
                TSeqsValue seqs_trim = trimPrefix(matches, lentokeep);
                SEQAN_OMP_PRAGMA(critical)
                {
                    assignValue(seqs, readPos, seqs_trim);
                }
            }
        }        
    },
    Standard(), threading);
}

// ----------------------------------------------------------------------------
// Function trim() - QualSld
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TIdents, typename TMap, typename TThreading>
inline void 
trim(TSeqs & seqs, TIdents & idents, TMap & map, TThreading const & threading,
QualSld)
{    
    typedef typename Iterator<TSeqs const, Standard>::Type      TSeqsIt;
    typedef typename Size<TSeqs>::Type                          TReadPos;
    typedef typename Value<TSeqs const>::Type                   TSeqsValue;
    typedef typename Size<TSeqsValue>::Type                     TSize;
    typedef typename Id<TSeqs>::Type                            TId;

    typedef typename TMap::mapped_type                          TParam;
    typedef float                                               TFValue;
    
    TParam windowsLength = map.at("windows_length");
    int meanQuality = static_cast<int>(map.at("mean_quality"));

    iterate(seqs, [&](TSeqsIt const & it)
    {
        TReadPos readPos = position(it, seqs);
        TId id = positionToId(seqs, readPos);
        TSeqsValue const & matches = value(it);
        TSize size = length(matches);

        // Discard record.
        if(size <  windowsLength)
        {
            idents.insert(id);
            return;
        }

        // Init var.
        TFValue currSum = 0.0f;
        TSize lentokeep = size;

        // Main loop.
        for(int i = size-1; i > -1 ; --i)
        {
            if((size - i) >= windowsLength)
            {
                if (i+windowsLength >= size)
                {
                    currSum += getQualityValue(matches[i]);
                }
                else
                {
                    currSum += getQualityValue(matches[i]) -
                    getQualityValue(matches[i + windowsLength]);
                }

                if((currSum / windowsLength) < meanQuality)
                {
                    lentokeep = i;
                }
            }
            else
            {
                currSum += getQualityValue(matches[i]);
            }
        }

        // Crop record while last base is minus than mean quality required.
        if(lentokeep != size)
        {
            bool isCrop = false;
            while(getQualityValue(matches[lentokeep]) < meanQuality 
                && lentokeep > 1)
            {
                --lentokeep;
                isCrop = true;
            }
            if (isCrop)
                ++lentokeep; // convert inclusive to exclusive
        }
        
        // Choose keeping or discarding record.
        if(lentokeep < 1 or lentokeep < windowsLength)
        {
            SEQAN_OMP_PRAGMA(critical)
            {
                idents.insert(id);
            }
        }
        else if(lentokeep<size)
        {
            TSeqsValue seqs_trim = trimPrefix(matches, lentokeep);
            SEQAN_OMP_PRAGMA(critical)
            {
                assignValue(seqs, readPos, seqs_trim);
            }
        }
    },
    Standard(), threading);    
}

// ----------------------------------------------------------------------------
// Function trim() - LenMin
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TIdents, typename TMap, typename TThreading>
inline void
trim(TSeqs & seqs, TIdents & idents, TMap & map, TThreading const & threading,
LenMin)
{
    typedef typename Iterator<TSeqs const, Standard>::Type      TSeqsIt;
    typedef typename Value<TSeqs const>::Type                   TSeqsValue;
    typedef typename Size<TSeqsValue>::Type                     TSize;
    typedef typename Id<TSeqs>::Type                            TId;
    typedef typename TMap::mapped_type                          TParam;

    TParam lenMin = map.at("len_min");

    iterate(seqs, [&](TSeqsIt const & it)
    {
        TSeqsValue const & matches = value(it);
        TSize size = length(matches);
        if (size <= lenMin)
        {
            TId id = positionToId(seqs, position(it, seqs));
            SEQAN_OMP_PRAGMA(critical)
            {
                idents.insert(id);
            }
        }
    },
    Standard(), threading);
}

// ----------------------------------------------------------------------------
// Function infoInit() - Info
// ----------------------------------------------------------------------------

template <typename TSize, typename TValue>
inline void infoInit(TSize & size, TValue & rest, TValue & steps)
{
    if (size <= InfoValues<>::getWindowSize())
    {
        rest = size;
    } 
    else
    {
        steps = ((size - InfoValues<>::getWindowSize()) / 
                InfoValues<>::getWindowStep()) + 1;
        rest = size - steps * InfoValues<>::getWindowStep();
        while (rest <= InfoValues<>::getWindowStep())
        {
            rest += InfoValues<>::getWindowStep();
            steps--;
        }
    }
}

// ----------------------------------------------------------------------------
// Function trim() - InfoDust
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TIdents, typename TMap, typename TThreading>
inline void 
trim(TSeqs & seqs, TIdents & idents, TMap & map, TThreading const & threading,
InfoDust)
{

    typedef typename Iterator<TSeqs const, Standard>::Type    TSeqsIt;
    typedef typename Value<TSeqs const>::Type    TSeqsValue;
    typedef typename Size<TSeqsValue>::Type        TSize;
    typedef typename Id<TSeqs>::Type        TId;

    typedef typename TMap::mapped_type TParam;
    typedef float TFValue;

    TParam cutoff = map.at("score");

    iterate(seqs, [&](TSeqsIt const & it)
    {
        TId id = positionToId(seqs, position(it, seqs));

        TSeqsValue const & matches = value(it);
        TSize size = length(matches);

        // Init vars.
        TParam rest = 0, steps = 0;
        TParam start = 0;
        TFValue score = 0.0 ;

        String<unsigned> kmerCounts;
        Accumulator<TFValue, AccuAverage> vals;

        infoInit(size, rest, steps);

        // Count Kmers.
        for(unsigned i = 0; i < steps; ++i)
        {
            start = i * InfoValues<>::getWindowStep();
            countKmers(kmerCounts, 
                        infix(matches, start, start+InfoValues<>::getWindowSize()),
                        3);
            score = 0.0;
            for(unsigned i = 0; i < length(kmerCounts); ++i)
            {
                if (kmerCounts[i] == 0)
                    continue;
                score += kmerCounts[i] * (kmerCounts[i] - 1) * 0.5;
            }
            push(vals, score * InfoValues<>::getByNum());
        }

        if(rest > 5)
        {
            start = steps * InfoValues<>::getWindowStep();
            countKmers(kmerCounts, infix(matches, start, start+rest), 3);
            score = 0.0;
            for(unsigned i = 0; i < length(kmerCounts); ++i)
            {
                if (kmerCounts[i] == 0)
                    continue;
                score += kmerCounts[i] * (kmerCounts[i] - 1) * 0.5;
            }
            push(vals, 
                ((score / (rest - 3)) *  (InfoValues<>::getWindowMax()/(rest-2))) );
        }
        else
        {
            push(vals, 31);
        }

        // Get mean.
        if(trunc(average(vals) * 100 / 31 ) > cutoff)
        {
            SEQAN_OMP_PRAGMA(critical)
            {
                idents.insert(id);
            }
        }

    },
    Standard(), threading);
}

// ----------------------------------------------------------------------------
// Function trim() - InfoN
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TIdents, typename TMap, typename TThreading>
inline void
trim(TSeqs & seqs, TIdents & idents, TMap & map, TThreading const & threading,
InfoN)
{
    typedef typename Iterator<TSeqs, Standard>::Type        TSeqsIt;
    typedef typename Id<TSeqs>::Type                        TId;
    typedef typename Value<TSeqs const>::Type               TSeqsValue;
    typedef typename Iterator<TSeqsValue, Standard>::Type   TSeqsValueIt;
    typedef typename TMap::mapped_type                      TParam;

    TParam score = map.at("score");

    iterate(seqs, [&](TSeqsIt const & it)
    {
        TId id = positionToId(seqs, position(it, seqs));
        TSeqsValue const & matches = value(it);

        // Init vars.
        TParam count = 0;

        // Main loop.
        for (TSeqsValueIt it = begin(matches); 
            it != end(matches, Standard()); 
            ++it)
        {
            if (ordValue(value(it)) == 4)
            {
                ++count;
                if (count >= score)
                {
                    break;
                }
            }

        }

        // Discard.
        if (count >= score)
        {
            SEQAN_OMP_PRAGMA(critical)
            {
                idents.insert(id);
            }
        }

    },
    Standard(), threading);
}

// ----------------------------------------------------------------------------
// Function trim()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TParam, typename TFormat_>
inline void trim(SeqStore<TSpec, TConfig> & me, TParam & params, Tag<TFormat_>, SequencingSingle)
{
    typedef Tag<TFormat_> TFormat;
    trim(me.seqs.i1, me.idents, params, 
        typename TConfig::TThreading(), 
        TFormat());
}

template <typename TSpec, typename TConfig, typename TParam, typename TFormat_>
inline void trim(SeqStore<TSpec, TConfig> & me, TParam & params, Tag<TFormat_>, SequencingPaired)
{
    typedef Tag<TFormat_> TFormat;
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {
            trim(me.seqs.i1, me.idents, params, 
            Serial(), 
            TFormat());
        }
        SEQAN_OMP_PRAGMA(section)
        {
            trim(me.seqs.i2, me.idents, params, 
            Serial(), 
            TFormat());
        }
    }
}

// --------------------------------------------------------------------------
// Function chooseTrimmer()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStringSet, 
    typename TFormat_>
inline void
chooseTrimmer(SeqStore<TSpec, TConfig> & me, TStringSet const & stringSet,
Tag<TFormat_>)
{
    typedef Tag<TFormat_> TFormat;

    for(auto & idtrimmer : stringSet)
    {
        if (std::strcmp(toCString(idtrimmer.name), 
            IdTrimmer<TFormat_>::VALUE[0]) == 0)
        {
            trim(me, idtrimmer.map, TFormat(), typename TConfig::TSequencing());
            return;
        }
    }
}

template <typename TSpec, typename TConfig, typename TStringSet, typename TTag>
inline void
chooseTrimmer(SeqStore<TSpec, TConfig> & me, TStringSet const &
stringSet, TagList<TTag, void> const)
{
    chooseTrimmer(me, stringSet, TTag());
}


template <typename TSpec, typename TConfig, typename TStringSet, typename TTag,
typename TSubList>
inline void
chooseTrimmer(SeqStore<TSpec, TConfig> & me, TStringSet &stringSet,
TagList<TTag, TSubList> const)
{
    chooseTrimmer(me, stringSet, TTag());
    chooseTrimmer(me, stringSet, TSubList());
}

template <typename TSpec, typename TConfig, typename TStringSet,
typename TTagList> 
inline void chooseTrimmer(SeqStore<TSpec, TConfig> & me, 
TStringSet & stringSet, TagSelector<TTagList> const)
{
    chooseTrimmer(me, stringSet, TTagList());
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void clear(Triple<TValue> & triple, TripleOne)
{
    clear(triple.i1);
}

template <typename TValue>
inline void clear(Triple<TValue> & triple, TripleTwo)
{
    clear(triple, TripleOne());
    clear(triple.i2);
}

template <typename TValue>
inline void clear(Triple<TValue> & triple, TripleThree)
{
    clear(triple, TripleTwo());
    clear(triple.i3);
}

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me, SeqStoreValue)
{
    clear(me.names, TripleThree());
    clear(me.seqs, TripleThree());
}

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me, SeqStoreRemove)
{
    clear(me.namesRemove, TripleThree());
    clear(me.seqsRemove, TripleThree());
}

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me, SeqStoreIdent)
{
    me.idents.clear();
}

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me)
{
    clear(me, SeqStoreValue());
    clear(me, SeqStoreRemove());
    clear(me, SeqStoreIdent());
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline bool empty(SeqStore<TSpec, TConfig> const & me)
{
    return empty(me.names.i1);
}

// ----------------------------------------------------------------------------
// Function size()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TValue>
inline void size(SeqStore<TSpec, TConfig> const & me, TValue & len)
{
    SEQAN_ASSERT_EQ(length(me.names.i1), length(me.names.i2));
    len += length(me.names.i1);
}

// ----------------------------------------------------------------------------
// Function statsDistributionReads()
// ----------------------------------------------------------------------------

template <typename TSeqs, typename TMap> inline void 
statsDistributionReads(TSeqs & seqs, TMap & map)
{
    typedef typename Iterator<TSeqs const, Standard>::Type      TSeqsIt;
    typedef typename Value<TSeqs const>::Type                   TSeqsValue;
    typedef typename Size<TSeqsValue>::Type                     TSize;

    iterate(seqs, [&](TSeqsIt const & it)
    {
        TSeqsValue const & matches = value(it);
        TSize size = length(matches);
        map[size] += 1;
    },
    Standard(), Serial());    
}

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileSpec,
typename TValue>
inline void
readRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Input, TFileSpec>> & fileIn,
TValue const & readBatch, FFastq, SequencingSingle)
{
    readRecords(me.names.i1, me.seqs.i1, fileIn.i1, readBatch);
}

template <typename TSpec, typename TConfig, typename TFileSpec,
typename TValue>
inline void 
readRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Input, TFileSpec>> & fileIn, 
TValue const & readBatch, FFastq, SequencingPaired)
{
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {
            readRecords(me.names.i1, me.seqs.i1, fileIn.i1, readBatch);
        }
        SEQAN_OMP_PRAGMA(section)
        {
            readRecords(me.names.i2, me.seqs.i2, fileIn.i2, readBatch);
        }
    }
}

template <typename TSpec, typename TConfig, typename TFileSpec,
typename TValue>
inline void
readRecords(SeqStore<TSpec, TConfig>,
Pair<FormattedFile<Fastq, Input, TFileSpec>>, TValue, FInterleaved,
SequencingSingle)
{
    return;
}

template <typename TSpec, typename TConfig, typename TFileSpec,
typename TValue>
inline void
readRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Input, TFileSpec>> & fileIn,
TValue const & readBatch, FInterleaved, SequencingPaired)
{
    typedef typename TConfig::TAlphabet                 TAlphabet;
    typedef typename TConfig::TAlphabetSequence         TAlphabetSequence;
    typedef typename TConfig::TSeqSpec                  TSeqSpec;

    typedef typename TConfig::TSeqsNameSpec             TSeqsNameSpec;
    typedef typename TConfig::TSeqsSpec                 TSeqsSpec;

    typedef StringSet<TAlphabet, TSeqsNameSpec>         TSeqName;
    typedef typename Iterator<TSeqName, Standard>::Type TSeqNameIt;
    typedef typename Value<TSeqName>::Type              TSeqNameValue;
    typedef typename Size<TSeqName>::Type               TSize;

    typedef String<TAlphabetSequence, TSeqSpec>         TSeq;
    typedef StringSet<TSeq, TSeqsSpec>                  TSeqs;
    typedef typename Iterator<TSeqs, Standard>::Type    TSeqsIt;

    // Read all.
    readRecords(me.names.i3, me.seqs.i3, fileIn.i1, readBatch);

    TSize lenRecords = length(me.names.i3) / 2;
    SEQAN_ASSERT_EQ(lenRecords%2, 0);

    // Remove specific interleaved extension name.
    String<CharString> needles;
    appendValue(needles, "\\1");
    appendValue(needles, "\\2");
    Pattern<String<CharString>, WuManber> pattern(needles);
    CharString srep("");
    iterate(me.names.i3, [&](TSeqNameIt const & it)
    {
        TSeqNameValue & name = value(it);
        Finder<CharString> finder(name);
        if (find(finder, pattern))
        {
            replace(name, beginPosition(finder), endPosition(finder), srep);
        }
    },
    Standard(), Serial());
    
    // Reserve.
    reserve(me.names.i1, lenRecords, Exact());
    reserve(me.seqs.i1, lenRecords, Exact());

    reserve(me.names.i2, lenRecords, Exact());
    reserve(me.seqs.i2, lenRecords, Exact());
        
    // Populate.
    TSeqNameIt itNames = begin(me.names.i3);
    TSeqsIt itSeqs = begin(me.seqs.i3);
    TSize size = length(me.names.i3);
    
    for(TSize step = 0; step < size; ++step, ++itNames, ++itSeqs)
    {
        
        if (step%2==0)
        {
            appendValue(me.names.i1, value(itNames));
            appendValue(me.seqs.i1, value(itSeqs));
        }
        else
        {
            appendValue(me.names.i2, value(itNames));
            appendValue(me.seqs.i2, value(itSeqs));
        }
    }
}

template <typename TSpec, typename TConfig, typename TFileSpec,
typename TValue>
inline void
readRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Input, TFileSpec>> & fileIn, TValue const& readBatch)
{
    readRecords(me, fileIn, readBatch, 
    typename TConfig::TInputFormat(), 
    typename TConfig::TSequencing());
}

// ----------------------------------------------------------------------------
// Function writeRecords()
// ----------------------------------------------------------------------------

template <typename TFileSpec, typename TNames, typename TSeqs>
inline void
writeRecords(TNames & names, TSeqs & seqs,
FormattedFile<Fastq, Output, TFileSpec> & fileOut, FInterleaved)
{
    typedef typename Value<TNames, 1>::Type                 TNamesValue;
    typedef typename Value<TNamesValue>::Type               TTNamesValue;
    typedef typename Size<TNamesValue>::Type                TSize;
    typedef typename Iterator<TNamesValue, Standard>::Type  TNamesIt;
    
    typedef typename Value<TSeqs, 1>::Type                  TSeqsValue;
    typedef typename Iterator<TSeqsValue, Standard>::Type   TSeqsIt;
    
    // Init values.
    TSize lenRecordsForward = length(names.i1);
    TSize lenRecordsReverse = length(names.i2);
    SEQAN_ASSERT_EQ(lenRecordsForward, lenRecordsReverse);
    TSize lenRecordsTotal = lenRecordsForward + lenRecordsReverse;

    // Specific interleaved.
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {  
            reserve(names.i3, lenRecordsTotal, Exact());
            TNamesIt itNamesForward = begin(names.i1, Standard());
            TNamesIt itNamesReverse = begin(names.i2, Standard());
        
            for (;
                itNamesForward != end(names.i1) || itNamesReverse != end(names.i2); 
                ++itNamesForward, ++itNamesReverse)
            {
                TTNamesValue & nameForward = value(itNamesForward);
                append(nameForward, "\\1");
                appendValue(names.i3, nameForward);

                TTNamesValue & nameReverse = value(itNamesReverse);
                append(nameReverse, "\\2");
                appendValue(names.i3, nameReverse);
            }
        }
        SEQAN_OMP_PRAGMA(section)
        {
 
            reserve(seqs.i3, lenRecordsTotal, Exact());
            TSeqsIt itSeqsForward = begin(seqs.i1, Standard());
            TSeqsIt itSeqsReverse = begin(seqs.i2, Standard());
            for (;
                itSeqsForward != end(seqs.i1) || itSeqsReverse != end(seqs.i2); 
                ++itSeqsForward, ++itSeqsReverse)
            {
                appendValue(seqs.i3, value(itSeqsForward));
                appendValue(seqs.i3, value(itSeqsReverse));
            }
        }
    }
    // Write.
    writeRecords(fileOut, names.i3, seqs.i3);
}


template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Output, TFileSpec>> & fileOut, FInterleaved)
{
    writeRecords(me.names, me.seqs, fileOut.i1, FInterleaved());
}


template <typename TFileSpec, typename TNames, typename TSeqs>
inline void
writeRecords(TNames const & names, TSeqs const & seqs, 
FormattedFile<Fastq, Output, TFileSpec> & fileOut, FFastq, SequencingSingle)
{
    writeRecords(fileOut, names.i1, seqs.i1);
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Output, TFileSpec>> & fileOut, FFastq, SequencingSingle)
{
    writeRecords(me.names, me.seqs, fileOut.i1, FFastq(), SequencingSingle());
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Output, TFileSpec>> & fileOut, FFastq, SequencingPaired)
{
    SEQAN_OMP_PRAGMA(parallel sections)
    {
        SEQAN_OMP_PRAGMA(section)
        {   
            writeRecords(fileOut.i1, me.names.i1, me.seqs.i1);
        }
        SEQAN_OMP_PRAGMA(section)
        {
            writeRecords(fileOut.i2, me.names.i2, me.seqs.i2);
        }
    }
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Output, TFileSpec>> & fileOut, FFastq)
{
    writeRecords(me, fileOut, FFastq(), typename TConfig::TSequencing());
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
Pair<FormattedFile<Fastq, Output, TFileSpec>> & fileOut)
{
    writeRecords(me, fileOut, typename TConfig::TOutputFormat());
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me,
FormattedFile<Fastq, Output, TFileSpec> & fileOut, SequencingPaired)
{
    writeRecords(me.namesRemove, me.seqsRemove, fileOut, FInterleaved());
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void 
writeRecords(SeqStore<TSpec, TConfig> & me,
FormattedFile<Fastq, Output, TFileSpec> & fileOut, SequencingSingle)
{
    writeRecords(me.namesRemove, me.seqsRemove, fileOut, FFastq(), 
    SequencingSingle());
}

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void
writeRecords(SeqStore<TSpec, TConfig> & me, 
FormattedFile<Fastq, Output, TFileSpec> & fileOut)
{
    writeRecords(me, fileOut, typename TConfig::TSequencing());
}

#endif  // #ifndef APP_HMNTRIMMER_SEQS_H_
