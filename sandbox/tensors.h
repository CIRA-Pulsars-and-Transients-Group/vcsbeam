#include "units.h"

/** Class for an axis */
template<class T>
class Axis
{
    private:
        std::vector<T> values;
        SIUnit u;

        T scale;
        T offset;
        long unsigned int size;

        bool use_custom_values;

    public:
        Axis();
}

/** Class for multi-dimensional tensors */
class Tensor
{
    private:
        void *data;   /** The pointer to the tensor data on the CPU */
        void *d_data; /** The pointer to the tensor data on the GPU */
}
