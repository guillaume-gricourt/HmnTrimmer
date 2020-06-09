// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Miscalleanous tags
// ============================================================================
#ifndef APP_HMNTRIMMER_MISCTAG_H_
#define APP_HMNTRIMMER_MISCTAG_H_

using namespace seqan;

// ============================================================================
// Enum
// ============================================================================

// ----------------------------------------------------------------------------
// Enum mode file
// ----------------------------------------------------------------------------

enum class FileStreamFormat : unsigned
{
    Undefined,
    Fastq,
    Interleaved
};

// ----------------------------------------------------------------------------
// Enum sequencing
// ----------------------------------------------------------------------------

enum class Sequencing : unsigned
{
    Undefined,
    Single,
    Paired
};

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tags format file
// ----------------------------------------------------------------------------

struct FFastq_;
typedef Tag<FFastq_>                FFastq;

struct FInterleaved_;
typedef Tag<FInterleaved_>          FInterleaved;

// ----------------------------------------------------------------------------
// Tags sequencing
// ----------------------------------------------------------------------------

struct SequencingSingle_;
typedef Tag<SequencingSingle_>      SequencingSingle;

struct SequencingPaired_;
typedef Tag<SequencingPaired_>      SequencingPaired;

// ----------------------------------------------------------------------------
// Tags SeqStore
// ----------------------------------------------------------------------------

struct SeqStoreIdent_;
typedef Tag<SeqStoreIdent_>         SeqStoreIdent;

struct SeqStoreRemove_;
typedef Tag<SeqStoreRemove_>        SeqStoreRemove;

struct SeqStoreValue_;
typedef Tag<SeqStoreValue_>         SeqStoreValue;

// ----------------------------------------------------------------------------
// Tags Triple
// ----------------------------------------------------------------------------

struct TripleOne_;
typedef Tag<TripleOne_>             TripleOne;

struct TripleTwo_;
typedef Tag<TripleTwo_>             TripleTwo;

struct TripleThree_;
typedef Tag<TripleThree_>           TripleThree;

#endif  // #ifndef APP_HMNTRIMMER_MISCTAG_H_
