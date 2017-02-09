/*!
 *  Allocates memory to hold data
 *  @tparam T template class
 *  @param[in] data container to hold data
 *  @param[in] n size of data
 */
template <class T>
void bitpit::vtk::allocate( T & data, int n){
    BITPIT_UNUSED(data);
    BITPIT_UNUSED(n);
}

/*!
 *  Allocates memory to hold data
 *  @tparam T template class
 *  @param[in] data container to hold data
 *  @param[in] n size of data
 */
template <class T>
void bitpit::vtk::allocate( std::vector<T> & data, int n){
    data.resize(n);
}
