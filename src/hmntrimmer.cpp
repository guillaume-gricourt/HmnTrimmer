// ============================================================================
//                                HmnTrimmer
// ============================================================================
//
// ============================================================================
// Author: Gricourt Guillaume guillaume.gricourt@aphp.fr
// ============================================================================
// Comment: Main file 
// ============================================================================
#define HMNTRIMMER

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// STL headers
// ----------------------------------------------------------------------------

#include <type_traits>
#include <random>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/arg_parse.h> 
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "hmntrimmer.hpp"
#include "version.hpp"

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void
setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "HmnTrimmer");
    setShortDescription(parser, "HmnTrimmer, a trimmer of NGS reads");
    setCategory(parser, "Trimmer");

    setDate(parser, HMNTRIMMER_DATE);
    setVersion(parser, HMNTRIMMER_VERSION);

    // Setup mandatory arguments.
    addUsageLine(parser, "HmnTrimmer [\\fIOPTIONS\\fP] [\\fITRIMMERS\\fP]");

    // Setup files.
    addSection(parser, "Input/Output options");

    //addText(parser, "Files Fastq standards - single/paired");
    addOption(parser, ArgParseOption("iff", "input-fastq-forward", "File with\
    read forward for apply trimmers.", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "input-fastq-forward", Options::getProcessExt());
    addOption(parser, ArgParseOption("ifr", "input-fastq-reverse", "File with\
    read reverse for apply trimmers.", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "input-fastq-reverse", Options::getProcessExt());
    addOption(parser, ArgParseOption("off", "output-fastq-forward", "File with\
    write forward trimmed.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-fastq-forward", Options::getProcessExt());
    addOption(parser, ArgParseOption("ofr", "output-fastq-reverse", "File with\
    write reverse trimmed.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-fastq-reverse", Options::getProcessExt());

    //addText(parser, "Files Fastq interleaved");
    addOption(parser, ArgParseOption("ifi", "input-fastq-interleaved", "File\
    with read interleaved", ArgParseOption::INPUT_FILE));
    setValidValues(parser, "input-fastq-interleaved", Options::getProcessExt());
    addOption(parser, ArgParseOption("ofi", "output-fastq-interleaved", "File\
    with write reverse trimmed.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-fastq-interleaved", Options::getProcessExt());    

    //addText(parser, "File to output discarded reads");
    addOption(parser, ArgParseOption("u", "output-fastq-discard", "File with\
    discard sequences.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-fastq-discard", Options::getProcessExt());
        
    // Setup trimmers.
    addSection(parser, "Trimmers. Classified in severals categories : quality,\
    length and information");
    //addText(parser, "Quality");
    addOption(parser, ArgParseOption("qualtail", "quality-tail", "",
    ArgParseOption::STRING));
    addOption(parser, ArgParseOption("slidingwindow", "quality-sliding-window",
    "", ArgParseOption::STRING));

    addOption(parser, ArgParseOption("lenmin", "length-min", "",
    ArgParseOption::INTEGER));

    addOption(parser, ArgParseOption("infodust", "information-dust", "",
    ArgParseOption::INTEGER));
    addOption(parser, ArgParseOption("infon", "information-n", "",
    ArgParseOption::INTEGER));

    // Setup performance options.
    addSection(parser, "Performance/Other Options");

    addOption(parser, ArgParseOption("r", "output-report", "File output report",
    ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output-report", Options::getReportExt());

    addOption(parser, ArgParseOption("t", "threads", "Specify the number of\
    threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
#ifdef _OPENMP
    setMaxValue(parser, "threads", "8");
#else
    setMaxValue(parser, "threads", "1");
#endif
    setDefaultValue(parser, "threads", options.threadsCount);

    addOption(parser, ArgParseOption("rb", "reads-batch", "Specify the \
    number of reads to process in one batch.", ArgParseOption::INTEGER));
    setMinValue(parser, "reads-batch", "100");
    setMaxValue(parser, "reads-batch", "50000000");
    setDefaultValue(parser, "reads-batch", options.readsBatch);

    addOption(parser, ArgParseOption("ver", "verbose", "Specify the log \
    level to use", ArgParseOption::INTEGER));
    setMinValue(parser, "verbose", "1");
    setMaxValue(parser, "verbose", "6");
    setDefaultValue(parser, "verbose", options.logLevel);
}

// ----------------------------------------------------------------------------
// Function parseStringArg()
// ----------------------------------------------------------------------------

void splitStringArg(CharString & sarg, String<CharString> & res)
{
    typedef Infix<CharString>::Type TInfix;
    typedef Suffix<CharString>::Type TSuffix;
    Finder<CharString> finder(sarg);
    Pattern<CharString, Horspool> p(":");

    unsigned begin = 0, end = 0;
    while (find(finder, p))
    {
        end = position(finder);
        TInfix inf(sarg, begin, end); 
        appendValue(res, CharString(inf));
        begin = end + 1;
    }
    if (length(res) > 0)
    {
        TSuffix suf(sarg, begin);
        appendValue(res, CharString(suf));
    }
    else
    {
        appendValue(res, sarg);
    }
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char
const ** argv)
{
    ArgTrimmer argTrimmer;
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Get files input.
    bool isInputFastqForward = isSet(parser, "input-fastq-forward");
    bool isInputFastqReverse = isSet(parser, "input-fastq-reverse");
    bool isInputFastqInterleaved = isSet(parser, "input-fastq-interleaved");
    if (isInputFastqForward and ! isInputFastqReverse and !
    isInputFastqInterleaved)
    {
        getOptionValue(options.inputFile.i1, parser, "input-fastq-forward");
        options.formatInput = FileStreamFormat::Fastq;
        options.sequencing = Sequencing::Single;
    }
    else if (isInputFastqForward and isInputFastqReverse and !
    isInputFastqInterleaved)
    {
        getOptionValue(options.inputFile.i1, parser, "input-fastq-forward");
        getOptionValue(options.inputFile.i2, parser, "input-fastq-reverse");
        options.formatInput = FileStreamFormat::Fastq;
        options.sequencing = Sequencing::Paired;
    }
    else if(! isInputFastqForward and ! isInputFastqReverse and
    isInputFastqInterleaved)
    {
        getOptionValue(options.inputFile.i1, parser, "input-fastq-interleaved");
        options.formatInput = FileStreamFormat::Interleaved;
        options.sequencing = Sequencing::Paired;
    }
    if (options.formatInput == FileStreamFormat::Undefined or options.sequencing
    == Sequencing::Undefined)
    {
        std::cerr << getAppName(parser) << "Incompatibles files in input are \
        indicated : either \"input-fastq-forward\", \"input-fastq-forward\" \
        and \"input-fastq-reverse\", \"input-fastq-interleaved\"" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Get files output.
    bool isOutputFastqForward = isSet(parser, "output-fastq-forward");
    bool isOutputFastqReverse = isSet(parser, "output-fastq-reverse");
    bool isOutputFastqInterleaved = isSet(parser, "output-fastq-interleaved");
    if (isOutputFastqForward and ! isOutputFastqReverse and !
    isOutputFastqInterleaved)
    {
        getOptionValue(options.outputFile.i1, parser, "output-fastq-forward");
        options.formatOutput = FileStreamFormat::Fastq;
    }
    else if (isOutputFastqForward and isOutputFastqReverse and !
    isOutputFastqInterleaved)
    {
        getOptionValue(options.outputFile.i1, parser, "output-fastq-forward");
        getOptionValue(options.outputFile.i2, parser, "output-fastq-reverse");
        options.formatOutput = FileStreamFormat::Fastq;
    }
    else if(! isOutputFastqForward and ! isOutputFastqReverse and
    isOutputFastqInterleaved)
    {
        getOptionValue(options.outputFile.i1, parser, "output-fastq-interleaved");
        options.formatOutput = FileStreamFormat::Interleaved;
    }
    if (options.formatOutput == FileStreamFormat::Undefined)
    {
        std::cerr << getAppName(parser) << "Incompatibles files in output \
        are indicated : either \"output-fastq-forward\", \
        \"output-fastq-forward\" and \"output-fastq-reverse\", \
        \"output-fastq-interleaved\"" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Check paired-end compatibility.
    bool isWellFormat = false;
    if (options.sequencing == Sequencing::Single and options.formatInput ==
    FileStreamFormat::Fastq and options.formatOutput == FileStreamFormat::Fastq)
        isWellFormat = true;
    if (options.sequencing == Sequencing::Paired and 
        (options.formatInput == FileStreamFormat::Fastq or options.formatInput
        == FileStreamFormat::Interleaved) and
        (options.formatOutput == FileStreamFormat::Fastq or options.formatOutput
        == FileStreamFormat::Interleaved))
        isWellFormat = true;
    if (!isWellFormat)
    {
        std::cerr << getAppName(parser) << "Conflict single/paired formats \
        between input/output files is found" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Ouput discarded records.
    if(getOptionValue(options.discardFile, parser, "output-fastq-discard"))
        options.isDiscardFile = true;

    // Output report.
    if(getOptionValue(options.reportFile, parser, "output-report"))
        options.isReportFile = true;

    // Get performance options.
    getOptionValue(options.threadsCount, parser, "threads");
    getOptionValue(options.readsBatch, parser, "reads-batch");
    getOptionValue(options.logLevel, parser, "verbose");

    // Get trimmers.
    CharString strimmer("");
    int itrimmer(0);
    
    if(getOptionValue(strimmer, parser, "quality-tail")){
        argTrimmer.setName(IdTrimmer<QualTail_>::VALUE[0]);
        splitStringArg(strimmer, argTrimmer.sarg);
        if(argTrimmer.sizeArg() == 1)
        {
            argTrimmer.setParam("base_quality", argTrimmer.sarg[0]);
            argTrimmer.setParam("base_number", CharString("1"));
        }
        else if(argTrimmer.sizeArg() == 2)
        {
            argTrimmer.setParam("base_quality", argTrimmer.sarg[0]);
            argTrimmer.setParam("base_number", argTrimmer.sarg[1]);
        }
        else if(argTrimmer.sizeArg() == 3)
        {
            argTrimmer.setParam("base_quality", argTrimmer.sarg[0]);
            argTrimmer.setParam("base_number", argTrimmer.sarg[1]);
            argTrimmer.setParam("len_perc", argTrimmer.sarg[2]);
        }
        else
        {
            std::cerr << getAppName(parser) << ": Quality Tail trimmer ->\
            required format \"<base quality>:<number of bases>\" (<int>:<int>)\
            or \"<base quality>:1\" (<int>)" << std::endl;
            return ArgumentParser::PARSE_ERROR;
        }
        appendValue(options.trimmers, argTrimmer);
        argTrimmer.clear();
    }
        
    if(getOptionValue(strimmer, parser, "quality-sliding-window")){
        argTrimmer.setName(IdTrimmer<QualSld_>::VALUE[0]);
        splitStringArg(strimmer, argTrimmer.sarg);

        if(argTrimmer.sizeArg() == 2)
        {
            argTrimmer.setParam("mean_quality", argTrimmer.sarg[0]);
            argTrimmer.setParam("windows_length", argTrimmer.sarg[1]);
        }
        else
        {
            std::cerr << getAppName(parser) << ": Quality Sliding Windows\
            trimmer -> required format \"<windows length>:<mean quality>\"\
            (<int>:<int>)" << std::endl;
            return ArgumentParser::PARSE_ERROR;
        }
        appendValue(options.trimmers, argTrimmer);
        argTrimmer.clear();
    }

    if(getOptionValue(itrimmer, parser, "length-min")){
        argTrimmer.setName(IdTrimmer<LenMin_>::VALUE[0]);
        argTrimmer.setParam("len_min", itrimmer);
        appendValue(options.trimmers, argTrimmer);
        argTrimmer.clear();
    }
    
    if(getOptionValue(itrimmer, parser, "information-dust")){
        argTrimmer.setName(IdTrimmer<InfoDust_>::VALUE[0]);
        argTrimmer.setParam("score", itrimmer);
        appendValue(options.trimmers, argTrimmer);
        argTrimmer.clear();
    }
    
    if(getOptionValue(itrimmer, parser, "information-n")){
        argTrimmer.setName(IdTrimmer<InfoN_>::VALUE[0]);
        argTrimmer.setParam("score", itrimmer);
        appendValue(options.trimmers, argTrimmer);
        argTrimmer.clear();
    }
    
    // Get software.
    options.softName = getAppName(parser);
    options.version = getVersion(parser);

    // Get command line.
    for (auto i = 0; i < argc; i++) {
        append(options.commandLine, argv[i]);
        appendValue(options.commandLine, ' ');
    }
    eraseBack(options.commandLine);

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function configureTrimming()
// ----------------------------------------------------------------------------

template <typename TThreading, typename TSequencing, typename TInputFormat>
void 
configureTrimming(Options const & options, TThreading const & threading, 
TSequencing const & sequencing, TInputFormat const & formatInput)
{
    switch (options.formatOutput)
    {
        case FileStreamFormat::Fastq:
            return spawnTrimming(options, threading, sequencing, formatInput, FFastq());
        case FileStreamFormat::Interleaved:
            return spawnTrimming(options, threading, sequencing, formatInput, FInterleaved());
        default:
            return;
    }
}

template <typename TThreading, typename TSequencing>
void 
configureTrimming(Options const & options, TThreading const & threading, 
TSequencing const & sequencing)
{
    switch (options.formatInput)
    {
        case FileStreamFormat::Fastq:
            return configureTrimming(options, threading, sequencing, FFastq());
        case FileStreamFormat::Interleaved:
            return configureTrimming(options, threading, sequencing,
            FInterleaved());
        default:
            return;
    }
}

template <typename TThreading>
void configureTrimming(Options const & options, TThreading const & threading)
{
    switch (options.sequencing)
    {
        case Sequencing::Single:
            return configureTrimming(options, threading, SequencingSingle());
        case Sequencing::Paired:
            return configureTrimming(options, threading, SequencingPaired());
        default:
            return;
    }
}

void configureTrimming(Options const & options)
{
#ifdef _OPENMP
    if (options.threadsCount > 1)
        configureTrimming(options, Parallel());
    else
#endif
        configureTrimming(options, Serial());
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    if (argc <= 1)
    {
        printHelp(parser, std::cout, "txt", true);
        return ArgumentParser::PARSE_ERROR;
    }

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc,
    argv);

    if (res != ArgumentParser::PARSE_OK)
    {
        printHelp(parser, std::cout, "txt", true);
        return ArgumentParser::PARSE_ERROR;
    }
    // Log.
    options.logger->info("Programme {}", toCString(getAppName(parser)));
    options.logger->info("Run : {}", toCString(options.commandLine));
    options.logger->debug("Using threads : {}", options.threadsCount);
    options.logger->debug("Log level : {}", options.logLevel);
    options.logger->debug("Reads batch : {}", options.readsBatch);

    try
    {
        configureTrimming(options);
    }
    catch (Exception const & e)
    {
        options.logger->error(e.what());
        return 1;
    }
    options.logger->info("Analyse is finished");
    return 0;
}
