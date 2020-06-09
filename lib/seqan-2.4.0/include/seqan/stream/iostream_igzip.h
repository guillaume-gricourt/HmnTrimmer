#ifndef INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_H_
#define INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_H_

#include <vector>
#include <iostream>
#include <algorithm>

namespace igzip_stream {

const size_t IGZIP_BLOCK_SIZE = 921600;//1048576; // 1024 * 1024
const size_t IGZIP_LEVEL_DEFAULT = 0;

// ===========================================================================
// Classes
// ===========================================================================

// --------------------------------------------------------------------------
// Class basic_igzip_streambuf
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_igzip_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &              ostream_reference;
    typedef ElemA                                       char_allocator_type;
    typedef ByteT                                       byte_type;
    typedef ByteAT                                      byte_allocator_type;
    typedef byte_type *                                 byte_buffer_type;
    typedef Tr                                          traits_type;
    typedef typename Tr::char_type                      char_type;
    typedef typename Tr::int_type                       int_type;
    typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
    typedef std::vector<char_type, char_allocator_type> char_vector_type;

    // Construct a igzip stream
    basic_igzip_streambuf(ostream_reference ostream_,
                        size_t level_,
                        size_t buffer_size_);

    ~basic_igzip_streambuf();

    int sync();
    int_type overflow(int_type c);

    // flushes the zip buffer and output buffer.
    // This method should be called at the end of the compression.
    // Calling flush multiple times, will lower the compression ratio.
    std::streamsize flush();

    // flushes the zip buffer and output buffer and finalize the zip stream
    // This method should be called at the end of the compression.
    std::streamsize flush_finalize();


private:
    bool igzip_to_stream(char_type *, std::streamsize);
    size_t fill_input_buffer();
    // flush the zip buffer using a particular mode and flush output buffer
    std::streamsize flush(int flush_mode);

    ostream_reference m_ostream;
    LZ_Stream2 m_igzip_stream;

    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_zip_ostreambase
// --------------------------------------------------------------------------
// Base class for igzip ostreams.
// Contains a basic_igzip_streambuf.

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_igzip_ostreambase :
    virtual public std::basic_ios<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &                      ostream_reference;
    typedef basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> igzip_streambuf_type;

    // Construct a zip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_igzip_ostreambase(ostream_reference ostream_,
                          size_t level_,
                          size_t buffer_size_) :
        m_buf(ostream_, level_, buffer_size_)
    {
        this->init(&m_buf);
    }

    // returns the underlying zip ostream object
    igzip_streambuf_type * rdbuf() { return &m_buf; }

private:
    igzip_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_zip_ostream
// --------------------------------------------------------------------------
// A igzipper ostream

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_igzip_ostream :
    public basic_igzip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
    public std::basic_ostream<Elem, Tr>
{
public:
    typedef basic_igzip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT> igzip_ostreambase_type;
    typedef std::basic_ostream<Elem, Tr>                          ostream_type;
    typedef ostream_type &                                        ostream_reference;

    // Constructs a igzipper ostream decorator
    //
    // ostream_ ostream where the compressed output is written
    // level_ level of compression
    // buffer_size_ level of compression

    basic_igzip_ostream(ostream_reference ostream_,
                      size_t level_ = IGZIP_LEVEL_DEFAULT,
                      size_t buffer_size_ = IGZIP_BLOCK_SIZE) :
        igzip_ostreambase_type(ostream_, level_, buffer_size_),
        ostream_type(this->rdbuf())
    {}

    ~basic_igzip_ostream()
    {
        //this->flush(); this->rdbuf()->flush_finalize();
    }

    // flush inner buffer and zipper buffer
    basic_igzip_ostream<Elem, Tr> & zflush()
    {
        this->flush(); this->rdbuf()->flush(); return *this;
    }
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_zip_ostream<char>
typedef basic_igzip_ostream<char>     igzip_ostream;
// A typedef for basic_zip_ostream<wchar_t>
typedef basic_igzip_ostream<wchar_t>  igzip_wostream;

} // namespace igzip_stream

#endif // INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_H_
