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

    Stats() :
        totalReads(0),
        keepReads(0),
        time(0)
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
// Function printStats()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void printStats(Trimming<TSpec, TConfig> const & me)
{
    // Init.
    rapidjson::Document document;
    document.SetObject();
    rapidjson::Document::AllocatorType& allocator = document.GetAllocator();

    // Fill.
    rapidjson::Value ksoftware(rapidjson::kObjectType);
    ksoftware.AddMember("name", 
        rapidjson::StringRef(toCString(me.options.softName)), 
        allocator);
    ksoftware.AddMember("version", 
        rapidjson::StringRef(toCString(me.options.version)), 
        allocator);
    document.AddMember("software", ksoftware, allocator);

    rapidjson::Value kstatistics(rapidjson::kObjectType);
    kstatistics.AddMember("runtime", me.stats.time, allocator);
    kstatistics.AddMember("total", me.stats.totalReads, allocator);
    kstatistics.AddMember("keep", me.stats.keepReads, allocator);
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
        size(me.reads, me.stats.totalReads);
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
        size(me.reads, me.stats.keepReads);
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
