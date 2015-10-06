// ========================================================================== //
//                         - SORTING ALGORITHMS -                             //
//                                                                            //
// Functions for data sorting.                                                //
// ========================================================================== //
// INFO                                                                       //
// Author    : Alessandro Alaia                                               //
// Version   : v2.0                                                           //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __SORT_ALG_HH__
# define __SORT_ALG_HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <iostream>

// Classes
// none

// Optimad C++ Library
# include "Operators.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //

// boolean vectors
typedef vector< bool >                 bvector1D;
typedef vector< bvector1D >            bvector2D;
typedef vector< bvector2D >            bvector3D;
typedef vector< bvector3D >            bvector4D;

// characters vectors
typedef vector< char >                 cvector1D;
typedef vector< cvector1D >            cvector2D;
typedef vector< cvector2D >            cvector3D;
typedef vector< cvector3D >            cvector4D;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// ========================================================================== //
// CLASSES                                                                    //
// ========================================================================== //

// kdtree ------------------------------------------------------------------- //

template <class T, class T1 = int>
class kdnode {

    // Members ========================================================== //
    public:
    int        lchild_;                                                   // pointer to left child
    int        rchild_;                                                   // pointer to left child
    T         *object_;                                                   // pointer to object
    T1         label;

    // Constructor ====================================================== //
    public:
    kdnode(                                                               // default constructor for kdnode variables
        void                                                              // (input) none
    );

    // Destructor ======================================================= //
    public:
    ~kdnode(                                                              // default destructor for kdnode variables
        void                                                              // (input) none
    );
};

template <int d, class T, class T1 = int>
class kdtree {

    // Members ============================================================== //
    public:
    int                     MAXSTK;                                           // max stack size
    int                     n_nodes;                                          // number of nodes
    vector< kdnode<T, T1> > nodes;                                            // kd-tree nodes

    // Constructors ========================================================= //
    public:
    kdtree(                                                                   // Default constructor for kdtree
        int                 stack_size = 10                                   // (input/optional stack size)
    );

    // Destructors ========================================================== //
    public:
    ~kdtree(                                                                  // default destructor for kdtree variables
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    int exist(                                                                // Check if element exist in the kd-tree
        T           *                                                         // (input) pointer to element to be tested
    );
    int exist(                                                                // Check if element exist in the kd-tree
        T           *,                                                        // (input) pointer to element to be tested
        T1          &                                                         // (input/output) label of the kd node matching test object
    );

    template <class T2>
    int h_neighbor(                                                           // Check if a kd-node exists in the h-neighborhood of a given item
        T           *,                                                        // (input) pointer to element to be tested
        T2           ,                                                        // (input) radius of ball
        bool         ,
        int         n = 0,                                                    // (input/optional) root for binary search algorithm
        int         l = 0                                                     // (input/optional) level of root on binary tree
    );

    template <class T2>
    int h_neighbor(                                                           // Check if a kd-node exists in the h-neighborhood of a given item
        T           *,                                                        // (input) pointer to element to be tested
        T1          &,                                                      // (input/output) label of the kd node matching test object
        T2           ,                                                        // (input) radius of ball
        int         n = 0,                                                    // (input/optional) root for binary search algorithm
        int         l = 0                                                     // (input/optional) level of root on binary tree
    );
    void insert(                                                              // Insert new element in the kd-tree
        T           *                                                         // (input) pointer to element to be inserted
    );
    void insert(                                                              // Insert new element in the kd-tree
        T           *,                                                        // (input) pointer to element to be inserted
        T1          &                                                         // (input) label of the new element
    );
    private:
    void IncreaseStack(                                                       // Increase stack size
        void                                                                  // )input) none
    );
    void DecreaseStack(                                                       // Decrease stack size
        void                                                                  // )input) none
    );

};

// min PQUEUE --------------------------------------------------------------- //
template <class T, class T1 = T>
class minPQUEUE {

    // Members ============================================================== //
    public:
    int                 MAXSTK;                                               // Maximal stack size between resize
    int                 heap_size;                                            // number of elements in stack
    vector< T >         keys;                                                 // stack
    vector< T1 >        labels;                                               // labels associated to keys
    ivector2D          *map;                                                  // pointer to mapper
    private:
    bool                use_labels;                                           // flag for key labelling

    // Constructor ========================================================== //
    public:
    minPQUEUE(                                                                // Default constructor for min priority queue
        bool            a = false,                                            // (input/optional) flag for key labelling
        ivector2D      *b = NULL                                              // (input/optional) pointer to user-defined map
    );
    minPQUEUE(                                                                // Default constructor for min priority queue
        int              ,                                                    // (input) stack size
        bool            a = false,                                            // (input/optional) flag for key labelling
        ivector2D      *b = NULL                                              // (input/optional) pointer to user-defined map
    );

    // Destructor =========================================================== //
    public:
    ~minPQUEUE(                                                               // Standard destructor for min priority queues
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear min heap content
        void                                                                  // (input) none
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &                                                     // (input/output) root value
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &,                                                    // (input/output) root value
        T1              &                                                     // (input/output) root label
    );
    void insert(                                                              // Insert a new key
        T               &                                                     // (input) new key to be inserted
    );
    void insert(                                                              // Insert a new key
        T               &,                                                    // (input) new key to be inserted
        T1              &                                                     // (input) label
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &                                                     // (input) new key value
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &,                                                    // (input) new key value
        T1              &                                                     // (input) label attached to the new key
    );
    void build_heap(                                                          // Build min-heap
        void                                                                  // (input) none
    );
    void display(                                                             // Display min-heap content
        ostream         &                                                     // (input) output stream
    );
    private:
    void IncreaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void DecreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (input) none
    );
    void heapify(                                                             // Restore min heap condition on spacified element
        int                                                                   // (input) position of element in stack
    );

};

// max PQUEUE --------------------------------------------------------------- //
template <class T, class T1 = T>
class maxPQUEUE {

    // Members ============================================================== //
    public:
    int                 MAXSTK;                                               // Maximal stack size between resize
    int                 heap_size;                                            // number of elements in stack
    vector< T >         keys;                                                 // stack
    vector< T1 >        labels;                                               // labels associated to keys
    ivector2D          *map;                                                  // pointer to mapper
    private:
    bool                use_labels;                                           // flag for key labelling

    // Constructor ========================================================== //
    public:
    maxPQUEUE(                                                                // Default constructor for min priority queue
        bool            a = false,                                            // (input/optional) flag for key labelling
        ivector2D      *b = NULL                                              // (input/optional) pointer to user-defined map
    );
    maxPQUEUE(                                                                // Default constructor for min priority queue
        int              ,                                                    // (input) stack size
        bool            a = false,                                            // (input/optional) flag for key labelling
        ivector2D      *b = NULL                                              // (input/optional) pointer to user-defined map
    );

    // Destructor =========================================================== //
    public:
    ~maxPQUEUE(                                                               // Standard destructor for min priority queues
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear max heap content
        void                                                                  // (input) none
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &                                                     // (input/output) root value
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &,                                                    // (input/output) root value
        T1              &                                                     // (input/output) root label
    );
    void insert(                                                              // Insert a new key
        T               &                                                     // (input) new key to be inserted
    );
    void insert(                                                              // Insert a new key
        T               &,                                                    // (input) new key to be inserted
        T1              &                                                     // (input) label
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &                                                     // (input) new key value
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &,                                                    // (input) new key value
        T1              &                                                     // (input) label attached to the new key
    );
    void build_heap(                                                          // Build min-heap
        void
    );
    void display(                                                             // Display min-heap content
        ostream         &                                                     // (input) output stream
    );
    private:
    void IncreaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void DecreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (input) none
    );
    void heapify(                                                             // Restore min heap condition on spacified element
        int                                                                   // (input) position of element in stack
    );

};

// LIFO stack --------------------------------------------------------------- //
template <class T>
class LIFOstack {

    // Members ============================================================== //
    public:
    int                 MAXSTK;                                               // Maximal stack size between resize
    int                 TOPSTK;                                               // Current stack size
    vector<T>           STACK;                                                // LIFO stack

    // Constructor ========================================================== //
    public:
    LIFOstack(                                                                // Standard constructor for LIFO stack
        void                                                                  // (input) none
    );
    LIFOstack(                                                                // Custom constructor #1 for LIFO stack
        int                                                                   // (input) maximal stack size
    );
    LIFOstack(                                                                // Custom constructor #2 for LIFO stack
        vector<T>       &                                                     // (input) items to be added in the LIFO stack
    );

    // Destructor =========================================================== //
    public:
    ~LIFOstack(                                                               // Standard destructor for LIFO stack
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear stack content
        void                                                                  // (input) none
    );
    void IncreaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void DecreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (output) none
    );
    T pop(                                                                    // Pop last item from stack
        void                                                                  // (input) none
    );
    void push(                                                                // Pusk item into the stack
        T                                                                     // (input) item to be pushed into the stack list
    );
    void push(                                                                // Pusk items into the stack
        vector<T>       &                                                     // (input) items to be pushed into the stack list
    );
    void display(                                                             // Display LIFO stack infos
        ostream         &                                                     // (input) output stream
    );
};

// FIFO stack --------------------------------------------------------------- //
template <class T>
class FIFOstack {

    // Members ============================================================== //
    public:
    int                 MAXSTK;                                               // Stack capacity
    int                 BEGSTK;                                               // Stack begin
    int                 ENDSTK;                                               // Stack end
    int                 TOPSTK;                                               // Current stack size
    vector<T>           STACK;                                                // FIFO stack

    // Constructor ========================================================== //
    public:
    FIFOstack(                                                                // Standard constructor for FIFO stack
        void                                                                  // (input) none
    );
    FIFOstack(                                                                // Custom constructor #1 for FIFO stack
        int                                                                   // (input) stack capacity
    );
    FIFOstack(                                                                // Custom constructor #2 for FIFO stack
        vector<T>       &                                                     // (input) items to be stocked in the FIFO stack
    );

    // Destructor =========================================================== //
    public:
    ~FIFOstack(                                                               // Standard destructor for LIFO stack
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void IncreaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void DecreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (input) none
    );
    int STACKempty(                                                           // Returns the number of empty position in the stack
        void                                                                  // (input) none
    );
    T pop(                                                                    // Pop last item from stack
        void                                                                  // (input) none
    );
    void push(                                                                // Push a single item into the stack
        T                                                                     // (input) item to be pushed into the stack list
    );
    void push(                                                                // Pusk multiple items into the stack
        vector<T>       &                                                     // (input) items to be pushed into the stack list
    );
    void clear(                                                               // Reset stack
        void                                                                  // (input) none
    );
    void display(                                                             // Display FIFO infos
        ostream         &                                                     // (input/output) output stream
    );
};

// ========================================================================== //
// FUNCTIONS PROTOTYPES                                                       //
// ========================================================================== //

// RANDOM SORTING =========================================================== //
void Extract_wo_Repl(                                                         // Extract integers without replacement
    int          ,                                                            // (input) number of integers to be extracted
    int          ,                                                            // (input) upper bound of extraction interval
    ivector1D   &                                                             // (input/output) list of extracted value
);

        // BIN SORTING ======================================================================= //
        template <class T>
        void Bin_Sort(                            // SORT DATA ON UNIFORM BINS
            int             ,                     // number of bins                            (INPUT)
            T               ,                     // data lower bound                          (INPUT)
            T               ,                     // data upper bound                          (INPUT)
            vector<T>       ,                     // data to be sorted                         (INPUT)
            ivector2D      &                      // list of data index for each bin           (INPUT/OUTPUT)
        );

        // HEAP SORT ========================================================================= //
        template <class T1, class T2>
        void max_Heapify(                         // RESTORE MAX-HEAP PROPERTY IN A BINARY TREE
            int             ,                     // Element to be restored                    (INPUT)
            int             ,                     // heap size                                 (INPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void min_Heapify(                         // RESTORE MAX-HEAP PROPERTY IN A BINARY TREE
            int             ,                     // Element to be restored                    (INPUT)
            int             ,                     // heap size                                 (INPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void Build_maxHeap(                       // BUILD MAX HEAP
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void Build_minHeap(                       // BUILD MAX HEAP
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void max_HeapSort(                        // SORT DATA IN DESCENDING ORDER USING A MAX-HEAP
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void min_HeapSort(                        // SORT DATA IN DESCENDING ORDER USING A MAX-HEAP
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *                      // mapper                                    (INPUT/OUTPUT)
        );

        // PRIORITY QUEUE ==================================================================== //
        template <class T1, class T2>
        void Extract_maxHeap(                     // EXTRACT ROOT FROM A MAX-HEAP TREE
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1             &,                     // extracted value                           (INPUT/OUTPUT)
            T2             &                      // data associated to the extracted value    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void Extract_minHeap(                     // EXTRACT ROOT FROM A MAX-HEAP TREE
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1             &,                     // extracted value                           (INPUT/OUTPUT)
            T2             &                      // data associated to the extracted value    (INPUT/OUTPUT)
        );
        template <class T1, class T2>
        void Increase_key(                        // INCREASE KEY VALUE IN A MAX-HEAP TREE
            int             ,                     // index of key to be increased              (INPUT)
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1              ,                     // new key value                             (OUTPUT)
            T2                                    // new data value attached to the key        (OUTPUT)
        );
        template <class T1, class T2>
        void Decrease_key(                        // DECREASE KEY VALUE IN A Min-HEAP TREE
            int             ,                     // index of key to be decreased              (INPUT)
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1              ,                     // new key value                             (OUTPUT)
            T2                                    // new data value attached to the key        (OUTPUT)
        );
        template <class T1, class T2>
        void Insert_maxHeap(                      // INSERT A NEW KEY IN A MAX-HEAP TREE
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1              ,                     // value of the new key                      (INPUT)
            T2                                    // data associated to the new key            (INPUT)
        );
        template <class T1, class T2>
        void Insert_minHeap(                      // INSERT A NEW KEY IN A MIN-HEAP TREE
            int            &,                     // heap size                                 (INPUT/OUTPUT)
            vector<T1>     &,                     // sorting keys                              (INPUT/OUTPUT)
            vector<T2>     *,                     // data attached to sorting keys             (INPUT/OUTPUT)
            ivector2D      *,                     // mapper                                    (INPUT/OUTPUT)
            T1              ,                     // value of the new key                      (INPUT)
            T2                                    // data associated to the new key            (INPUT)
        );

        // =================================================================================== //
        // TEMPLATES                                                                           //
        // =================================================================================== //
        # include "LIFOstack.tpp"
        # include "FIFOstack.tpp"
        # include "PQUEUE.tpp"
        # include "kdtree.tpp"

        # endif
