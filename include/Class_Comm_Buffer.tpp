/*
 * Class_Comm_Buff.tpp
 *
 *  Created on: 19/mar/2014
 *      Author: Marco Cisternino
 */


template<class T>
void Class_Comm_Buffer::write(T& val) {
	MPI_Datatype datatype = convert<T>();
	int error = MPI_Pack(&val,1,datatype,commBuffer,commBufferSize,&pos,MPI_COMM_WORLD);
}

template<class T>
void Class_Comm_Buffer::read(T& val) {
	MPI_Datatype datatype = convert<T>();
	int error = MPI_Unpack(commBuffer,commBufferSize,&pos,&val,1,datatype,MPI_COMM_WORLD);
}

