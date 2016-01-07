#ifndef __PATCHMAN_BINARY_STREAM_HPP__
#define __PATCHMAN_BINARY_STREAM_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <stdexcept>

// Bitpit
// none

// Others
// none

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //
// none

// Forward declarations ------------------------------------------------- //
class ibinarystream;
class obinarystream;

// Function prototypes -------------------------------------------------- //
template<typename T>
ibinarystream& operator>>(                                                  // Stream operator for class ibinarystream
    ibinarystream                     &istm,                                // (input) input stream
    T                               &val                                  // (input) value to be streamed
);
template<>
ibinarystream& operator>>(                                                  // Explicit specialization of input stream operator for std::string
    ibinarystream                     &istm,                                // (input) input stream
    std::string                     &val                                  // (input) string to be streamed
);
template<typename T>
obinarystream& operator<<(                                                  // Stream operator for class obinarystream
    obinarystream                     &ostm,                                // (input) output stream
    const T                         &val                                  // (input) value to be streamed
);
template<>
obinarystream& operator<<(                                                  // Explicit specialization of input stream operator for std::string
    obinarystream                     &ostm,                                // (input) output stream
    const std::string               &val                                  // (input) string to be streamed
);
obinarystream& operator<<(                                                  // Stream operator for class obinarystream
    obinarystream                     &ostm,                                // (input) output stream
    const char                      *val                                  // (input) char array to be streamed
);

// Class ibinarystream ---------------------------------------------------- //
class ibinarystream {

    // Member(s) ======================================================== //
    private:

    std::vector<char>               buffer;                               // stream buffer
    size_t                          current_pos;                          // Cursor position

    // Constructor(s) =================================================== //
    public:

    ibinarystream(                                                          // Default constructor (empty stream)
        void                                                              // (input) none
    );
    ibinarystream(                                                          // Custom constructor #1 (empty stream with known size);
        size_t                       size                                 // (input) buffer size
    );
    ibinarystream(                                                          // Custom constructor #2 (stream pointing to memory location)
        const char*                  buf_,                                // (input) pointer to memory location
        size_t                       size                                 // (input) buffer size
    );
    ibinarystream(                                                          // Custom constructor #3 (stream initialized from std::vector<char>)
        const std::vector<char>          &vec                                  // (input) vector used for initialization
    );

    // Destructor(s) ==================================================== //
    // default

    // Assignament operator(s) ========================================== //
    const ibinarystream& operator=(
        const ibinarystream           &istm
    ) = delete;

    // Public method(s) ================================================= //
    public:
    void open(                                                            // Open input stream from memory location
        const char                  *mem,                                 // (input) pointer to memory location
        size_t                       size                                 // (input) size (in bytes) of memory chunk
    );
    void close(                                                           // Close input stream from memory
        void                                                              // (input) none
    );
    bool eof(                                                             // Flag for eof
        void                                                              // (input) none
    ) const;
    std::ifstream::pos_type tellg(                                             // Returns current position of cursor in the buffer
        void                                                              // (input) none
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        size_t                       pos                                  // (input) position
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        std::streamoff               offset,                              // (input) offset with respect to the specified direction
        std::ios_base::seekdir       way                                  // (input) offset direction
    );
    const std::vector<char>& get_internal_vec(                                 // Returns reference to buffer
        void                                                              // (input) none
    ) { return(buffer); }
    char* get_buffer(                                                     // Returns pointer to buffer
        void                                                              // (input) none
    ) { return( buffer.data() ); }

    // Private methods(s) =============================================== //
    private:

    template<typename T>
    void read(                                                            // Read data from memory location pointed by t and store into stream buffer
        T                           &t                                    // (input) data to be imported in the stream buffer
    );
    void read(                                                            // Explicit template specialization of ibinarystream::read for vector<char>
        std::vector<char>           &vec                                  // (input) source vector
    );
    void read(                                                                // Read data from memory location pointed by p and store into stream buffer
        char                        *p,                                   // (input) pointer to memory location
        size_t                       size                                 // (input) size (in bytes) of data to be read
    );

    // Friendships ====================================================== //
    template< typename T >
    friend ibinarystream& operator >> (ibinarystream&, T& );
    friend ibinarystream& operator >> (ibinarystream&, std::string& );
};

// Class obinarystream ---------------------------------------------------- //
class obinarystream {

    // Member(s) ======================================================== //
    private:

    size_t                           current_pos;                         // Cursor current position
    std::vector<char>                buffer;                              // Buffer

    // Constructor(s) =================================================== //
    public:

    obinarystream(                                                          // Default constructor (create empty object)
        void                                                              // (input) none
    );
    obinarystream(                                                          // Custom constructor #1 (create an empty object with buffer of specified size)
        size_t                       size                                 // (input) none
    );

    // Destructor(s) ==================================================== //
    // none

    // Assignement operator(s) ========================================== //
    const obinarystream& operator=(
        const obinarystream           &
    ) = delete;

    // Public method(s) ================================================= //
    void open(                                                            // Open output stream
        size_t                       size                                 // (input) stream size
    );
    void close(                                                           // Close output stream
        void                                                              // (input) none
    );
    bool eof(                                                             // Flag for eof
        void                                                              // (input) none
    ) const;
    std::ifstream::pos_type tellg(                                             // Returns current position of cursor in the buffer
        void                                                              // (input) none
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        size_t                       pos                                  // (input) position
    );
    bool seekg (                                                          // Set cursor position in the current buffer
        std::streamoff               offset,                              // (input) offset with respect to the specified direction
        std::ios_base::seekdir       way                                  // (input) offset direction
    );
    const std::vector<char>& get_internal_vec(                                 // Returns reference to buffer
        void                                                              // (input) none
    ) { return(buffer); }
    char* get_buffer(                                                     // Returns pointer to buffer
        void                                                              // (input) none
    ) { return( buffer.data() ); }

    // Private method(s) ================================================ //
    private:

    template<typename T>
    void write(                                                           // Write data to internal buffer
        const T                     &t                                    // (input) data to be written in the internal buffer
    );
    void write(                                                           // Write char array to internal buffer
        const char                  *p,                                   // (input) pointer to char array
        size_t                       size                                 // (input) size of data chunk to be written in the internal buffer
    );
    void write(                                                           // Write vector of char to internal buffer
        const std::vector<char>     &vec                                  // (input) vector of char to be written in the internal buffer
    );

    // Friendship(s) ==================================================== //
    template<typename T>
    friend obinarystream& operator<<( obinarystream&, const T& );
    friend obinarystream& operator<<( obinarystream&, const std::string& );
    friend obinarystream& operator<<( obinarystream&, const char* );
};

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "binary_stream.tpp"

#endif
