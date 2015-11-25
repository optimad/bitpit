
// Operator "sign" ------------------------------------------------------------------- //
template <class T>
T sign(                                                                               // RETURNS THE SIGN 
        const  T                & val                                                      // (input) input value
      ){
    return (T(0) < val) - (val < T(0)) ;
};



