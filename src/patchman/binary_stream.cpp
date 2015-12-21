// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include "binary_stream.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// IMPLEMENTATIONS OF METHODS FOR CLASS ibinarystream                         //
// ========================================================================== //

// Constructor(s) =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        Default constructor. Initialize an empty object of class ibinarystream

*/
ibinarystream::ibinarystream(
    void
) {
    current_pos = 0;
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #1. Initialize an empty object of class ibinarystream
        with assigned size

        \param[in] size buffer size
*/
ibinarystream::ibinarystream(
    size_t                      size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    buffer.resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #2. Initialize a object of class ibinarystream
        pointing to a memory location with specified size

        \param[in] buf_ pointer to memory location
        \param[in] size buffer size

*/
ibinarystream::ibinarystream(
    const char                  *buf_,
    size_t                       size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    buffer.assign(buf_, buf_ + size);
}

// -------------------------------------------------------------------------- //
/*!
        Custom constructor #3. Initialize a object of class ibinarystream
        from a std::vector of char

        \param[in] vec input vector

*/
ibinarystream::ibinarystream(
    const vector<char>          &vec
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(vec.size());
    buffer.assign(vec.begin(), vec.end());
}

// Destructor(s) ============================================================ //
// default

// Assignament operator(s) ================================================== //
// disabled

// Public methods =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        Open stream from memory

        \param[in] mem pointer to memory location
        \param[in] size size (in bytes) of memory location to be streamed

*/
void ibinarystream::open(
    const char                  *mem,
    size_t                       size
) {
    current_pos = 0;
    buffer.clear();
    buffer.reserve(size);
    buffer.assign(mem, mem + size);
}

// -------------------------------------------------------------------------- //
/*!
        Close stream from memory

*/
void ibinarystream::close(
    void
)
{
    buffer.clear();
}

// -------------------------------------------------------------------------- //
/*!
        Returns true if end of file condition is met.

        \result boolean flag (true) if end of file is reached, (false) otherwise

*/
bool ibinarystream::eof(
    void
) const
{
    return current_pos >= buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Returns cursor position within the current buffer

        \result cursor position

*/
ifstream::pos_type ibinarystream::tellg(
    void
) {
    return current_pos;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] pos cursor new position

        \result (true) if new position is valid, (false) otherwise

*/
bool ibinarystream::seekg (
    size_t                       pos
) {
    if(pos<buffer.size())
        current_pos = pos;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] offset new offset position
        \param[in] way offset direction

        \result (true) if new position is valid, (false) otherwise

*/
bool ibinarystream::seekg (
    streamoff                    offset,
    ios_base::seekdir            way
) {
    if ( ( way == ios_base::beg ) && ( offset < (long) buffer.size() ) )
        current_pos = offset;
    else if ( ( way == ios_base::cur ) && ( current_pos + offset < buffer.size() ) )
        current_pos += offset;
    else if ( ( way == ios_base::end ) && ( (long) buffer.size() - offset >= 0 ) )
        current_pos = buffer.size() - offset;
    else
        return false;

    return true;
}

// Private method(s) ======================================================== //

// -------------------------------------------------------------------------- //
/*!
        Read data from memory location pointed by p and store into stream buffer

        \param[in] p pointer to memory location where data are stored
        \param[in] size size (in bytes) of data chunk to be read

*/
void ibinarystream::read(
    char                        *p,
    size_t                       size
) {
    if ( eof() || (current_pos + size) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>( p ), &buffer[current_pos], size);
    current_pos += size;
}

// -------------------------------------------------------------------------- //
/*!
        Read data from std::vector<char> and store content into stream buffer

        \param[in] vec vector of char

*/
void ibinarystream::read(
    std::vector<char>           &vec
) {
    if ( eof() || ( current_pos + vec.size() ) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>(&vec[0]), &buffer[current_pos], vec.size());
    current_pos += vec.size();
}

// ========================================================================== //
// IMPLEMENTATIONS OF METHODS FOR CLASS obinarystream                           //
// ========================================================================== //

// Constructor(s) =========================================================== //

// -------------------------------------------------------------------------- //
/*!
        Default constructor. Initialize an empty object
*/
obinarystream::obinarystream(
    void
) {
    current_pos = 0;
}

// -------------------------------------------------------------------------- //
/*!
        Default constructor. Initialize an empty object with buffer of specified
        size.

        \param[in] size buffer size
*/
obinarystream::obinarystream(
    size_t                       size
) {
    current_pos = 0;
    open(size);
}

// Destructor(s) ============================================================ //
// default

// Assignament operator(s) ================================================== //
// disabled

// Public method(s) ========================================================= //

// -------------------------------------------------------------------------- //
/*!
        Open output stream

        \param[in] size stream size

*/
void obinarystream::open(
    size_t                       size
) {
    buffer.reserve(size);
    buffer.resize(size);
}

// -------------------------------------------------------------------------- //
/*!
        Close output stream

*/
void obinarystream::close(
    void
) {
    buffer.clear();
}

// -------------------------------------------------------------------------- //
/*!
        Returns true if end of file condition is met.

        \result boolean flag (true) if end of file is reached, (false) otherwise

*/
bool obinarystream::eof(
    void
) const
{
    return current_pos >= buffer.size();
}

// -------------------------------------------------------------------------- //
/*!
        Returns cursor position within the current buffer

        \result cursor position

*/
ifstream::pos_type obinarystream::tellg(
    void
) {
    return current_pos;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] pos cursor new position

        \result (true) if new position is valid, (false) otherwise

*/
bool obinarystream::seekg (
    size_t                       pos
) {
    if(pos < buffer.size())
        current_pos = pos;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Set cursor position within the current buffer

        \param[in] offset new offset position
        \param[in] way offset direction

        \result (true) if new position is valid, (false) otherwise

*/
bool obinarystream::seekg (
    streamoff                    offset,
    ios_base::seekdir            way
) {
    if ( ( way == ios_base::beg ) && ( offset < (long) buffer.size() ) )
        current_pos = offset;
    else if ( ( way == ios_base::cur ) && ( current_pos + offset < buffer.size() ) )
        current_pos += offset;
    else if ( ( way == ios_base::end ) && ( (long) buffer.size() - offset >= 0 ) )
        current_pos = buffer.size() - offset;
    else
        return false;

    return true;
}

// -------------------------------------------------------------------------- //
/*!
        Write char array to internal buffer

        \param[in] p pointer to char array
        \param[in] size size of memory chunck to be writting to the internal
        buffer

*/
void obinarystream::write(
    const char                  *p,
    size_t                       size
) {
    if ( buffer.size() - current_pos < size ) {
        buffer.resize( size + current_pos );
    }
    for (size_t i = 0;  i < size; ++i) {
        buffer[current_pos] = p[i];
        ++current_pos;
    } // next i
}

// -------------------------------------------------------------------------- //
/*!
        Write vector of char to internal buffer

        \param[in] vec vector of char to be written in the internal buffer
        buffer

*/

void obinarystream::write(
    const vector<char>          &vec
) {
    if ( buffer.size() - current_pos < vec.size() ) {
        buffer.resize( vec.size() + current_pos );
    }
    for(size_t i = 0; i < vec.size(); ++i) {
        buffer[current_pos] = vec[i];
        ++current_pos;
    } //next i
}

// ========================================================================== //
// OPERATORS
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Explicit template specialization for stream operator for class ibinarystream

        \param[in] istm input stream
        \param[in] val std::string to be streamed

*/
template<>
ibinarystream& operator>>(
        ibinarystream             &istm,
        string                  &val)
{
    int                 size = 0;

    istm.read(size);

    if(size<=0)         return istm;

    std::vector<char> vec((size_t)size);
    istm.read(vec);
    val.assign(&vec[0], (size_t)size);

    return istm;
}

// -------------------------------------------------------------------------- //
/*!
        Stream std::string to internal buffer.

        \param[in] ostm output stream
        \param[in] val string

*/
template<>
obinarystream& operator << (
    obinarystream                 &ostm,
    const string                &val
) {
    int size = val.size();

    ostm.write(size);

    if(val.size()<=0)
        return ostm;

    ostm.write(val.c_str(), val.size());

    return ostm;
}

// -------------------------------------------------------------------------- //
/*!
        Stream char array to internal buffer.

        \param[in] ostm output stream
        \param[in] val pointer to char array

*/
obinarystream& operator<<(
    obinarystream                 &ostm,
    const char                  *val
) {
    int size = strlen(val);

    ostm.write(size);

    if(size<=0)
        return ostm;

    ostm.write(val, size);

    return ostm;
}
