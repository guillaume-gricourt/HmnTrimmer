#ifndef INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_IMPL_H_
#define INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_IMPL_H_

#include "iostream_igzip.h"

namespace igzip_stream {

/*
PLAN :
  - constructeur : ok
  - destructeur : ~
  - sync : ok
  - overflow :
  - igzip_to_stream :
  - flush() : 
  - flush(int) : 
*/

// --------------------------------------------------------------------------
// Class basic_zip_streambuf
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_igzip_streambuf(
    ostream_reference ostream_,
    size_t ,
    size_t buffer_size_
    ) :
    m_ostream(ostream_),
    m_output_buffer(buffer_size_, 0),
    m_buffer(buffer_size_, 0)
{
    init_stream(&m_igzip_stream);
    m_igzip_stream.end_of_stream = 0;

    this->setp(&(m_buffer[0]), &(m_buffer[m_buffer.size() - 1]));
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_igzip_streambuf()
{
    flush();
    m_ostream.flush();
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
int basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::sync()
{
    if (this->pptr() && this->pptr() > this->pbase())
    {
        if (traits_type::eq_int_type(overflow(traits_type::eof()), traits_type::eof()))
            return -1;
    }
    return 0;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
typename basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::overflow(
    typename basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type c)
{
    int w = static_cast<int>(this->pptr() - this->pbase());
    if (!traits_type::eq_int_type(c, traits_type::eof()))
    {
        *this->pptr() = c;
        ++w;
    }

    if (igzip_to_stream(this->pbase(), w))
    {
        this->setp(this->pbase(), this->epptr() - 1);
        return c;
    }
    else
    {
        return traits_type::eof();
    }
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
bool basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::igzip_to_stream(
    typename basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::char_type * buffer_,
    std::streamsize buffer_size_)
{
    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_igzip_stream.next_in = (byte_buffer_type)buffer_;
    m_igzip_stream.avail_in = static_cast<uInt>(buffer_size_ * sizeof(char_type));
    m_igzip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_igzip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    do
    {
        fast_lz(&m_igzip_stream);

        written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_igzip_stream.avail_out;
        total_written_byte_size += written_byte_size;
        // ouput buffer is full, dumping to ostream
        m_ostream.write((const char_type *) &(m_output_buffer[0]),
                        static_cast<std::streamsize>(written_byte_size / sizeof(char_type)));

        // checking if some bytes were not written.
        if ((remainder = written_byte_size % sizeof(char_type)) != 0)
        {
            // copy to the beginning of the stream
            std::memmove(&(m_output_buffer[0]),
                &(m_output_buffer[written_byte_size - remainder]),
                remainder);
        }
        m_igzip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
        m_igzip_stream.next_out = &m_output_buffer[remainder];
    }
    while (m_igzip_stream.avail_out == 0);

    return true;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_igzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush()//int flush_mode)
{
    int const buffer_size = static_cast<int>(this->pptr() - this->pbase()); // amount of data currently in buffer

    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_igzip_stream.next_in = (byte_buffer_type) this->pbase();
    m_igzip_stream.avail_in = static_cast<uInt>(buffer_size * sizeof(char_type));
    m_igzip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_igzip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    m_igzip_stream.end_of_stream = 1;

    do
    {
        fast_lz(&m_igzip_stream);

        written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_igzip_stream.avail_out;
        total_written_byte_size += written_byte_size;
        // ouput buffer is full, dumping to ostream
        m_ostream.write((const char_type *) &(m_output_buffer[0]),
                        static_cast<std::streamsize>(written_byte_size / sizeof(char_type)));

        // checking if some bytes were not written.
        if ((remainder = written_byte_size % sizeof(char_type)) != 0)
        {
            // copy to the beginning of the stream
            std::memmove(&(m_output_buffer[0]),
                &(m_output_buffer[written_byte_size - remainder]),
                remainder);
        }
        m_igzip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
        m_igzip_stream.next_out = &m_output_buffer[remainder];
    }
    while (m_igzip_stream.avail_out == 0);

    m_ostream.flush();

    return total_written_byte_size;
}

} // namespace igzip_stream

#endif // INCLUDE_SEQAN_STREAM_IOSTREAM_IGZIP_IMPL_H_
