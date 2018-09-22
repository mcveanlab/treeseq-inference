#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include "bgen.h"

typedef struct {
    PyObject_HEAD
    struct bgen_file *bgen;
    struct bgen_vi *bgen_index;
    size_t num_variants;
    size_t num_samples;
} BgenReader;

static void
BgenReader_dealloc(BgenReader* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
    if (self->bgen != NULL) {
        bgen_close(self->bgen);
        self->bgen = NULL;
    }
    if (self->bgen_index != NULL) {
        bgen_free_index(self->bgen_index);
        self->bgen_index = NULL;
    }
}

static int
BgenReader_init(BgenReader *self, PyObject *args, PyObject *kwds)
{
    int ret = 0;
    static char *kwlist[] = {"filename", NULL};
    struct bgen_var *variants = NULL;
    char *filename;

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &filename)) {
        ret = -1;
        goto out;
    }
    self->bgen = bgen_open(filename);
    if (self->bgen == NULL) {
        PyErr_SetString(PyExc_ValueError, "error opening file");
        ret = -1;
        goto out;
    }
    self->num_variants = bgen_nvariants(self->bgen);
    self->num_samples = bgen_nsamples(self->bgen);

    /* Only doing this to initialise the index. Not sure if this is the correct approach
     * to be honest. */
    variants = bgen_read_variants_metadata(self->bgen, &self->bgen_index, 0);
    if (variants == NULL) {
        PyErr_SetString(PyExc_ValueError, "error reading index");
        goto out;
    }
out:
    if (variants != NULL) {
        bgen_free_variants_metadata(self->bgen, variants);
    }
    return ret;
}

static PyObject *
BgenReader_get_probabilities(BgenReader* self, PyObject *args)
{
    PyObject *ret = NULL;
    unsigned long variant;
    struct bgen_vg *vg = NULL;
    npy_intp dims[2];
    PyArrayObject *probabilities = NULL;
     
    if (!PyArg_ParseTuple(args, "k", &variant)) {
        goto out;
    }
    if (variant >= self->num_variants) {
        PyErr_SetString(PyExc_ValueError, "Variant index out of bounds.");
        goto out;
    }
    vg = bgen_open_variant_genotype(self->bgen_index, variant);
    if (vg == NULL) {
        PyErr_SetString(PyExc_ValueError, "Error getting variant.");
        goto out;
    }
    dims[0] = self->num_samples;
    dims[1] = bgen_ncombs(vg);
    probabilities = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (probabilities == NULL) {
        goto out;
    }
    bgen_read_variant_genotype(self->bgen_index, vg, PyArray_DATA(probabilities));
    ret = (PyObject*) probabilities;
    probabilities = NULL;
out:
    if (vg != NULL) {
        bgen_close_variant_genotype(self->bgen_index, vg);
    }
    Py_XDECREF(probabilities);
    return ret;
}

static PyObject *
BgenReader_get_num_samples(BgenReader *self, void *closure)
{
    return Py_BuildValue("n", (Py_ssize_t) self->num_samples);
}

static PyObject *
BgenReader_get_num_variants(BgenReader *self, void *closure)
{
    return Py_BuildValue("n", (Py_ssize_t) self->num_variants);
}

static PyMethodDef BgenReader_methods[] = {
    {"get_probabilities", (PyCFunction)BgenReader_get_probabilities, METH_VARARGS,
        "Return the probabilities for the specified variant."},
    {NULL}  /* Sentinel */
};

static PyGetSetDef BgenReader_getsetters[] = {
    {"num_samples", (getter)BgenReader_get_num_samples, NULL, "number of samples", NULL},
    {"num_variants", (getter)BgenReader_get_num_variants, NULL, "number of variants", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject BgenReaderType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "simplebgen.BgenReader",             /* tp_name */
    sizeof(BgenReader), /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)BgenReader_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "BgenReader objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    BgenReader_methods,             /* tp_methods */
    0,             /* tp_members */
    BgenReader_getsetters,     /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)BgenReader_init,      /* tp_init */
    0,                         /* tp_alloc */
    0,                         /* tp_new */
};

static PyModuleDef simplebgenmodule = {
    PyModuleDef_HEAD_INIT,
    "simplebgen",
    "Example module that creates an extension type.",
    -1,
    NULL, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_simplebgen(void)
{
    PyObject* m;

    BgenReaderType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&BgenReaderType) < 0) {
        return NULL;
    }
    import_array();
    m = PyModule_Create(&simplebgenmodule);
    if (m == NULL) {
        return NULL;
    }
    Py_INCREF(&BgenReaderType);
    PyModule_AddObject(m, "BgenReader", (PyObject *)&BgenReaderType);
    return m;
}
