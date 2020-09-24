// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Main header 
// ============================================================================
#ifndef APP_HMNTRIMMER_H_
#define APP_HMNTRIMMER_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

// ----------------------------------------------------------------------------
// rapidjson headers
// ----------------------------------------------------------------------------

#include <rapidjson/document.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>

// ----------------------------------------------------------------------------
// spdlog headers
// ----------------------------------------------------------------------------

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "misc_tags.hpp"
#include "seqs.hpp"
#include "timer.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class RunConfig
// ----------------------------------------------------------------------------

template <typename TThreading_          = Parallel,
            typename TSequencing_       = SequencingSingle,
            typename TInputFormat_      = Fastq,
            typename TOutputFormat_     = Fastq>
struct RunConfig
{
    typedef TThreading_                 TThreading;
    typedef TSequencing_                TSequencing;
    typedef TInputFormat_               TInputFormat;
    typedef TOutputFormat_              TOutputFormat;
};

// ----------------------------------------------------------------------------
// Class ArgTrimmer
// ----------------------------------------------------------------------------

struct ArgTrimmer
{
    typedef CharString                  TName;
    typedef std::string                 TPName;
    typedef unsigned                    TValue;
    typedef std::map<TPName, TValue>    TMap; 
    typedef String<TName>               TArg;
    
    TName       name;
    TMap        map;
    TArg        sarg;

    ArgTrimmer() : 
        name(""),
        map(),
        sarg()
    {};

    void setName(TName const & sname)
    {
        name = sname;
    }

    bool setParam(TName const & sname, TName const & svalue)
    {
        TValue value(0);
        TPName name = std::string(toCString(sname));
        try
        {
            value = std::stoi(toCString(svalue));
        }
        catch (const std::exception& e)
        {
            return false;
        }
        map[name] = value;
        return true;
    }

    bool setParam(TName const & sname, int const & i)
    {
        TValue value(0);
        TPName name = std::string(toCString(sname));
        try
        {
            value = (unsigned) i;
        }
        catch (const std::exception& e)
        {
            return false;
        }
        map[name] = value;
        return true;
    }

    std::string arg2string() const
    {
        std::stringstream ss;
        for(auto const& arg : map)
        {
            ss << arg.first << ":" << arg.second << ", ";
        }
        std::string arg = ss.str();
        arg = arg.substr(0, arg.size()-2);
        return arg;
    }

    void clear ()
    {
        name = "";
        map.clear();
        seqan::clear(sarg);
    }

    unsigned sizeArg () const
    {
        return length(sarg);
    }
};

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{ 
    Pair<CharString>        inputFile;
    Pair<CharString>        outputFile;
    CharString              discardFile;
    CharString              reportFile;

    FileStreamFormat        formatInput;
    FileStreamFormat        formatOutput;
    Sequencing              sequencing;

    bool                    isDiscardFile;
    bool                    isReportFile;

    long                    readsBatch;
    unsigned                threadsCount;

    CharString              commandLine;
    CharString              version;
    CharString              softName;

    String<ArgTrimmer>      trimmers;

    std::shared_ptr<spdlog::logger>        logger;
    unsigned                logLevel;

    Options() :
        inputFile(),
        outputFile(),
        discardFile(""),
        reportFile(""),
        formatInput(FileStreamFormat::Undefined),
        formatOutput(FileStreamFormat::Undefined),
        sequencing(Sequencing::Undefined),
        isDiscardFile(false),
        isReportFile(false),
        readsBatch(1000000),
        threadsCount(1),
        logLevel(4)
    {
        logger = spdlog::stdout_logger_st("console");
        logger->set_pattern("%d-%m-%Y %R - %l - %v");
        switch (logLevel) {
            case 1: logger->set_level(spdlog::level::critical); break;
            case 2: logger->set_level(spdlog::level::err); break;
            case 3: logger->set_level(spdlog::level::warn); break;
            case 4: logger->set_level(spdlog::level::info); break;
            case 5: logger->set_level(spdlog::level::debug); break;
            case 6: logger->set_level(spdlog::level::trace); break;
            default: logger->set_level(spdlog::level::info); break;
        }
    }

    static constexpr const char* getProcessExt()
    {
        return "fq fastq fq.gz fastq.gz";
    }
    static constexpr const char* getReportExt()
    {
        return "json";
    }
};

// ----------------------------------------------------------------------------
// Class TrimmingTraits
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct TrimmingTraits
{
    //Parameters
    typedef typename TConfig::TThreading        TThreading;
    typedef typename TConfig::TSequencing       TSequencing;
    typedef typename TConfig::TInputFormat      TInputFormat;
    typedef typename TConfig::TOutputFormat     TOutputFormat;

    //Files
    typedef Pair<SeqFileIn>                     TReadsFileIn;
    typedef Pair<SeqFileOut>                    TReadsFileOut;
    typedef SeqFileOut                          TReadsFileDiscard;
};

// ----------------------------------------------------------------------------
// Class Stats
// ----------------------------------------------------------------------------

template <typename TValue>
struct Stats
{
    TValue  totalReads;
    TValue  keepReads;

    TValue  time;

    std::map<TValue, TValue>    distriBefore;
    std::map<TValue, TValue>    distriAfter;

    Stats() :
        totalReads(0),
        keepReads(0),
        time(0),
        distriBefore(),
        distriAfter()
    {}
};

// ----------------------------------------------------------------------------
// Class Trimming
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig = void>
struct Trimming
{
    typedef TrimmingTraits<TSpec, TConfig>          Traits;

    Options const &                                 options;
    Stats<unsigned long>                            stats;
    Timer<>                                         timer;

    //Store reads
    typedef SeqStore<void, SeqConfig<TConfig>>      TReads;
    TReads                                          reads;

    typename Traits::TReadsFileIn                   readsFileIn;
    typename Traits::TReadsFileOut                  readsFileOut;
    typename Traits::TReadsFileDiscard              readsFileDiscard;

    Trimming(Options const & options) :
        options(options)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function configureThreads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void configureThreads(Trimming<TSpec, TConfig> & me)
{
    omp_set_num_threads(me.options.threadsCount);
}

// ----------------------------------------------------------------------------
// Function printStatsMap()
// ----------------------------------------------------------------------------

template <typename TDocument, typename TJsonValue, typename TMap>
inline void printStatsMap(TJsonValue & jsonValue, TMap & map, TDocument & document)
{
    typedef typename TMap::mapped_type  TParam;
    typedef typename std::map<TParam, TParam>::const_iterator TIt;
    
    for (TIt it = map.begin(); it!=map.end(); ++it)
    {
        std::string skey = std::to_string(it->first);
        rapidjson::Value key(skey.c_str(), document.GetAllocator());
        rapidjson::Value value(static_cast<uint64_t>(it->second));
        jsonValue.AddMember(key, value, document.GetAllocator());
    }
}
 
// ----------------------------------------------------------------------------
// Function printStats()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void printStats(Trimming<TSpec, TConfig> const & me)
{
    typedef CharString                  TString;
    
    // Init.
    rapidjson::Document document;
    document.SetObject();
    rapidjson::Document::AllocatorType& allocator = document.GetAllocator();

    // Runtime.
    rapidjson::Value ksoftware(rapidjson::kObjectType);
    ksoftware.AddMember("name", 
        rapidjson::StringRef(toCString(me.options.softName)), 
        allocator);
    ksoftware.AddMember("version", 
        rapidjson::StringRef(toCString(me.options.version)), 
        allocator);
    document.AddMember("software", ksoftware, allocator);

    // Analyze.
    rapidjson::Value kanalyse(rapidjson::kObjectType);

    rapidjson::Value kruntime(rapidjson::kObjectType);
    kruntime.AddMember("unit", "seconds", allocator);
    kruntime.AddMember("value", me.stats.time, allocator);
    kanalyse.AddMember("runtime", kruntime, allocator);
 
    TString ssequencing = "";
    switch(me.options.sequencing){
        case Sequencing::Undefined : ssequencing = "undefined"; break;
        case Sequencing::Single : ssequencing = "single"; break;
        case Sequencing::Paired : ssequencing = "paired"; break;
    }
    kanalyse.AddMember("sequencing", 
        rapidjson::StringRef(toCString(ssequencing)), 
        allocator);
    
    rapidjson::Value kfilenames(rapidjson::kObjectType);
    rapidjson::Value kfilenamesInput(rapidjson::kArrayType);
    kfilenamesInput.PushBack(
        rapidjson::Value(toCString(me.options.inputFile.i1), allocator).Move(), 
        allocator);
    kfilenamesInput.PushBack(
        rapidjson::Value(toCString(me.options.inputFile.i2), allocator).Move(), 
        allocator);    
    kfilenames.AddMember("input", kfilenamesInput, allocator);

    rapidjson::Value kfilenamesOutput(rapidjson::kArrayType);
    kfilenamesOutput.PushBack(
        rapidjson::Value(toCString(me.options.outputFile.i1), allocator).Move(), 
        allocator);
    kfilenamesOutput.PushBack(
        rapidjson::Value(toCString(me.options.outputFile.i2), allocator).Move(), 
        allocator);    
    kfilenames.AddMember("output", kfilenamesOutput, allocator);
    kanalyse.AddMember("file", kfilenames, allocator);
    
    rapidjson::Value ktrimmers(rapidjson::kObjectType);
    for(auto & idtrimmer : me.options.trimmers)
    {
        ktrimmers.AddMember(
            rapidjson::StringRef(toCString(idtrimmer.name)),
            rapidjson::Value(idtrimmer.arg2string().c_str(), allocator),
            allocator);
    }
    kanalyse.AddMember("trimmers", ktrimmers, allocator);
    
    document.AddMember("analyze", kanalyse, allocator);

    // Statistics general.
    rapidjson::Value kstatistics(rapidjson::kObjectType);
    kstatistics.AddMember("total", me.stats.totalReads, allocator);
    kstatistics.AddMember("kept", me.stats.keepReads, allocator);
    kstatistics.AddMember("discarded", 
    me.stats.totalReads - me.stats.keepReads, allocator);

    rapidjson::Value kdistributionBefore(rapidjson::kObjectType);
    printStatsMap(kdistributionBefore, me.stats.distriBefore, document);
    kstatistics.AddMember("length_reads_before", 
    kdistributionBefore, allocator);

    rapidjson::Value kdistributionAfter(rapidjson::kObjectType);
    printStatsMap(kdistributionAfter, me.stats.distriAfter, document);
    kstatistics.AddMember("length_reads_after", 
    kdistributionAfter, allocator);
    
    document.AddMember("statistics", kstatistics, allocator);

    // Write output.
    std::ofstream ofs(toCString(me.options.reportFile));
    rapidjson::OStreamWrapper ofw(ofs); 
    rapidjson::Writer<rapidjson::OStreamWrapper> writer(ofw);
    document.Accept(writer);
}

// ----------------------------------------------------------------------------
// Function open() - open paired file
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool open(Pair<FormattedFile<TFileType, TDirection, TSpec> > & me,
        const char * fileName1,
        const char * fileName2)
{
    return open(me.i1, fileName1) && open(me.i2, fileName2);
}

// ----------------------------------------------------------------------------
// Function openInputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void
_openReadsInput(Trimming<TSpec, TConfig> & me, SequencingSingle, FFastq)
{
    if (!open(me.readsFileIn.i1, toCString(me.options.inputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsInput(Trimming<TSpec, TConfig> & me, SequencingSingle, FInterleaved)
{
    if (!open(me.readsFileIn.i1, toCString(me.options.inputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsInput(Trimming<TSpec, TConfig> & me, SequencingPaired, FFastq)
{
    if (!open(me.readsFileIn, 
        toCString(me.options.inputFile.i1), 
        toCString(me.options.inputFile.i2)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsInput(Trimming<TSpec, TConfig> & me, SequencingPaired, FInterleaved)
{
    if (!open(me.readsFileIn.i1, toCString(me.options.inputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void openInputFile(Trimming<TSpec, TConfig> & me)
{
    _openReadsInput(me, 
        typename TConfig::TSequencing(),
        typename TConfig::TInputFormat());
}

// ----------------------------------------------------------------------------
// Function openOutputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void
_openReadsOutput(Trimming<TSpec, TConfig> & me, SequencingSingle, FFastq)
{
    if (!open(me.readsFileOut.i1, toCString(me.options.outputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsOutput(Trimming<TSpec, TConfig> & me, SequencingSingle, FInterleaved)
{
    if (!open(me.readsFileOut.i1, toCString(me.options.outputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsOutput(Trimming<TSpec, TConfig> & me, SequencingPaired, FFastq)
{
    if (!open(me.readsFileOut, 
        toCString(me.options.outputFile.i1), 
        toCString(me.options.outputFile.i2)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void
_openReadsOutput(Trimming<TSpec, TConfig> & me, SequencingPaired, FInterleaved)
{
    if (!open(me.readsFileOut.i1, toCString(me.options.outputFile.i1)))
        throw RuntimeError("Error while opening reads file.");
}

template <typename TSpec, typename TConfig>
inline void openOutputFile(Trimming<TSpec, TConfig> & me)
{
    _openReadsOutput(me, 
        typename TConfig::TSequencing(), 
        typename TConfig::TOutputFormat());
}

// ----------------------------------------------------------------------------
// Function openDiscardFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void openDiscardFile(Trimming<TSpec, TConfig> & me)
{
    if (!open(me.readsFileDiscard, toCString(me.options.discardFile)))
        throw RuntimeError("Error while opening reads file.");
}

// ----------------------------------------------------------------------------
// Function close() - close single file
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline void 
close(Pair<FormattedFile<TFileType, TDirection, TSpec> > & me, 
SequencingSingle, FFastq)
{
    close(me.i1);
}

template <typename TFileType, typename TDirection, typename TSpec>
inline void 
close(Pair<FormattedFile<TFileType, TDirection, TSpec> > & me, 
SequencingSingle, FInterleaved)
{
    close(me.i1);
}

template <typename TFileType, typename TDirection, typename TSpec>
inline void
close(Pair<FormattedFile<TFileType, TDirection, TSpec> > & me, 
SequencingPaired, FInterleaved)
{
    close(me.i1);
}

// ----------------------------------------------------------------------------
// Function close() - close paired file
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline void 
close(Pair<FormattedFile<TFileType, TDirection, TSpec> > & me, 
SequencingPaired, FFastq)
{
    close(me.i1);
    close(me.i2);
}

// ----------------------------------------------------------------------------
// Function closeInputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void closeInputFile(Trimming<TSpec, TConfig> & me)
{
    close(me.readsFileIn, 
    typename TConfig::TSequencing(), 
    typename TConfig::TInputFormat());
}

// ----------------------------------------------------------------------------
// Function closeOutputFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void closeOutputFile(Trimming<TSpec, TConfig> & me)
{
    close(me.readsFileOut, 
    typename TConfig::TSequencing(), 
    typename TConfig::TOutputFormat());
}

// ----------------------------------------------------------------------------
// Function closeDiscardFile()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void closeDiscardFile(Trimming<TSpec, TConfig> & me)
{
    close(me.readsFileDiscard);
}

// ----------------------------------------------------------------------------
// Function loadReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadReads(Trimming<TSpec, TConfig> & me)
{
    readRecords(me.reads, me.readsFileIn, me.options.readsBatch);
}

// ----------------------------------------------------------------------------
// Function writeReadsDiscard()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void writeReadsDiscard(Trimming<TSpec, TConfig> & me)
{
    writeRecords(me.reads, me.readsFileDiscard);
}

// ----------------------------------------------------------------------------
// Function statsDistribution()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void statsDistributionReads(Trimming<TSpec, TConfig> & me, Before)
{
    statsDistributionReads(me.reads.seqs.i1, me.stats.distriBefore);
    if(IsSameType<typename TConfig::TSequencing, SequencingPaired>::VALUE)
    {
        statsDistributionReads(me.reads.seqs.i2, me.stats.distriBefore);
    }
}

template <typename TSpec, typename TConfig>
inline void statsDistributionReads(Trimming<TSpec, TConfig> & me, After)
{
    statsDistributionReads(me.reads.seqs.i1, me.stats.distriAfter);
    if(IsSameType<typename TConfig::TSequencing, SequencingPaired>::VALUE)
    {
        statsDistributionReads(me.reads.seqs.i2, me.stats.distriAfter);
    }
}

// ----------------------------------------------------------------------------
// Function writeReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void writeReads(Trimming<TSpec, TConfig> & me)
{
    writeRecords(me.reads, me.readsFileOut);
}

// ----------------------------------------------------------------------------
// Function chooseTrimmer()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void trim(Trimming<TSpec, TConfig> & me)
{
    chooseTrimmer(me.reads, me.options.trimmers, Trimmers());
}

// ----------------------------------------------------------------------------
// Function clearReads()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clearReads(Trimming<TSpec, TConfig> & me)
{
    clear(me.reads);
}

// ----------------------------------------------------------------------------
// Function runTrimming()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runTrimming(Trimming<TSpec, TConfig> & me)
{
    me.options.logger->trace("Start timer");
    start(me.timer);

    // Configure threads.
    me.options.logger->trace("Configure threading with OpenMP");
    configureThreads(me);

    // Open.
    me.options.logger->info("Open files Input");
    openInputFile(me);
    me.options.logger->trace("Open files Output");
    openOutputFile(me);
    if (me.options.isDiscardFile)
    {
        openDiscardFile(me);
    }
    // Process reads in blocks.
    unsigned batch = 1;
    while (true)
    {
        me.options.logger->debug("Batch : {}", batch);
        me.options.logger->debug("\tLoad");
        loadReads(me);
        if (me.options.isReportFile)
        {
            size(me.reads, me.stats.totalReads);
            statsDistributionReads(me, Before());
        }
        if (empty(me.reads)) break;
        me.options.logger->debug("\tTrim");
        trim(me);
        if (me.options.isDiscardFile)
        {
            me.options.logger->debug("\tUpdate discard");
            updateDiscard(me.reads);
            me.options.logger->debug("\tWrite discard");
            writeReadsDiscard(me);
        }
        me.options.logger->debug("\tUpdate");
        update(me.reads);
        if (me.options.isReportFile)
        {
            size(me.reads, me.stats.keepReads);
            statsDistributionReads(me, After());
        }
        me.options.logger->debug("\tWrite");
        writeReads(me);
        me.options.logger->debug("\tClear");
        clear(me.reads);

        // Update.
        ++batch;
    }

    // Close files.
    me.options.logger->debug("\tClose files");
    closeInputFile(me);
    closeOutputFile(me);
    me.options.logger->debug("\tClose discard files");
    if (me.options.isDiscardFile)
    {
        closeDiscardFile(me);
    }
    me.options.logger->trace("Stop timer");
    stop(me.timer);
    getTimer(me.timer, me.stats.time);

    // Write report.
    if(me.options.isReportFile)
    {
        me.options.logger->info("Write report : {}", 
        toCString(me.options.reportFile));
        printStats(me);
    }
}

// ----------------------------------------------------------------------------
// Function spawnTrimming()
// ----------------------------------------------------------------------------

template <typename TThreading, typename TSequencing, 
typename TInputFormat, typename TOutputFormat>
inline void spawnTrimming(Options const & options, TThreading const &,
TSequencing const &, TInputFormat const &, TOutputFormat const &)
{
    typedef RunConfig<TThreading, TSequencing, 
    TInputFormat, TOutputFormat>  TConfig;

    Trimming<void, TConfig> trimming(options);
    runTrimming(trimming);
}

#endif  // #ifndef APP_HMNTRIMMER_H_
