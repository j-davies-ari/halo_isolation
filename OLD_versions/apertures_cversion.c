#include <Python.h>
#include <stdbool.h>

#include <numpy/ndarraytypes.h>
#include <numpy/ndarrayobject.h>


void getPair(double *a, double *b) {    // Check that these are not pointing to NULL    assert(a);    assert(b);    *a = 1;    *b = getString();}



bool isvalueNOTinarray(int val, int *arr, int size){    int i;    for (i=0; i < size; i++) {        if (arr[i] == val)            return false;    }    return true;}

void c_aperture(double *xs, double *ys, double *zs, double *centre, double *sidesize, 
                double *boxsize, int n){ 
    
    sidesize = sidesize / *2.0;
    
    if (centre[0]+sidesize > boxsize ||
        centre[1]+sidesize > boxsize ||
        centre[2]+sidesize > boxsize)
        {

        xs = xs - (centre[0] - boxsize/2.0)
        xs = xs % boxsize
        xs = xs - boxsize/2.0

        ys = ys - (centre[1] - boxsize/2.0)
        ys = ys % boxsize
        ys = ys - boxsize/2.0

        zs = zs - (centre[2] - boxsize/2.0)
        zs = zs % boxsize
        zs = zs - boxsize/2.0

        }
    else
        {
        xs = xs - center[0]
        ys = ys - center[1]
        zs = zs - center[2]
        }

    int m;
    int mask [n];
    for(m=0;m<n;m++){ // Set these to values that don't make sense so that the wasted bit of the array can be chopped off later
        mask[m] = -1
        }

    int i;
    int m_ind;
    for(i=0;i<n;i++){
        m_ind = 0;
        if ((xs[i]**2+ys[i]**2+zs[i]**2)**0.5 < sidesize){
            mask[m_ind] = i;
            m_ind = m_ind + 1;
            }
        }
    
    for(i=0;i<n;i++){ // If a particle isn't in the mask, set its position to a big no. of Mpc for disposal
        if (isvalueNOTinarray(i,mask,n)){
            xs[i] = 100.;
            ys[i] = 100.;
            zs[i] = 100.;
            }
    }


static PyObject *aperture_module(PyObject *self, PyObject *args){
    // The inputs from python    PyArrayObject *xs_obj, *ys_obj, *zs_obj; 
    PyArrayObject *centre_obj;
    double *sidesize, *boxsize;

    // The C variables
    struct *xs, *ys, *zs;
    double *centre;
    int n;    if(!PyArg_ParseTuple(args, "OOOOii",&xs_obj, &ys_obj, &zs_obj, &centre_obj, &sidesize, &boxsize))
        return NULL;

    // point the C arrays to the numpy arrays
    xs = (double *)xs_obj->data;
    ys = (double *)ys_obj->data;
    zs = (double *)zs_obj->data;
    centre = (double *)centre_obj->data;

    // How many particles are there?
    n = (int) xs->dimensions[0];

    int mask [n];
    // do the work
    c_aperture(xs,ys,zs,centre,sidesize,boxsize,n)
    return Py_BuildValue("(o&,o&,o&,o&)",xs,ys,zs,mask);
}




static PyMethodDef ApertureMethods[] = {
  {"render", aperture_module, METH_VARARGS, "Method for creating apertures"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initaperture(void) {
  (void) Py_InitModule("aperture", ApertureMethods);
}
