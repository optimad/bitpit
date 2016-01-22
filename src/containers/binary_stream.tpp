// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS ibinarystream                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        \ingroup    BinaryStream
        Read data from memory location and store in the stream buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void ibinarystream::read(
    T                           &t
) {
    if ( eof() || ( current_pos + sizeof(T) ) > buffer.size() ) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>( &t ), &buffer[current_pos], sizeof(T));
    current_pos += sizeof(T);
}

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR CLASS obinarystream                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
/*!
        \ingroup    BinaryStream
        Write data from memory to internal buffer.

        \param[in] t data to be stored in the stream buffer

*/
template<typename T>
void obinarystream::write(
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
        \ingroup    BinaryStream
        Input stream operator for class ibinarystream.

        \param[in] istm input stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
ibinarystream& operator>> (
    ibinarystream                 &istm,
    T                           &val
) {
    istm.read(val);
    return istm;
}

// -------------------------------------------------------------------------- //
/*!
        \ingroup    BinaryStream
        Output stream operator for class ibinarystream.

        \param[in] ostm output stream
        \param[in] val  value to be streamed

        \result updated input stream

*/
template<typename T>
obinarystream& operator<< (
    obinarystream                 &ostm,
    const T                     &val
) {
    ostm.write(val);

    return ostm;
}

