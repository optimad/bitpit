// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS IBinaryStream                             //
// ========================================================================== //

namespace bitpit{

// -------------------------------------------------------------------------- //
/*!
        Read data from memory location and store in the stream buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void IBinaryStream::read(
    T                           &t
) {
    if ( eof() || ( current_pos + sizeof(T) ) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>( &t ), &buffer[current_pos], sizeof(T));
    current_pos += sizeof(T);
}

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS OBinaryStream                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Write data from memory to internal buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void OBinaryStream::write(
    const T                     &t
) {
    std::vector<char> vec(sizeof(T));
    memcpy( reinterpret_cast<void*>( &vec[0] ), reinterpret_cast<const void*>( &t ), sizeof(T) );
    write(vec);
}

// ========================================================================== //
// OPERATORS                                                                  //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        Input stream operator for class IBinaryStream.

        \param[in] istm input stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
IBinaryStream& operator>> (
    IBinaryStream                 &istm,
    T                           &val
) {
    istm.read(val);
    return istm;
}

// -------------------------------------------------------------------------- //
/*!
        Output stream operator for class IBinaryStream.

        \param[in] ostm output stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
OBinaryStream& operator<< (
    OBinaryStream                 &ostm,
    const T                     &val
) {
    ostm.write(val);

    return ostm;
}

}
