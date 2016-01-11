
// Operator "sign" ------------------------------------------------------------------- //
/*!
    Sign function.
    Given a a variable of integral type, val, returns:
    1 if val > 0
    -1, othersize.

    Template parameters can be any integral type such that operator< is defined.

    \result returns the sign of the input value.
*/
template <class T>
T sign(                                                                               // RETURNS THE SIGN 
        const  T                & val                                                      // (input) input value
      ){
    return (T(0) < val) - (val < T(0)) ;
};



