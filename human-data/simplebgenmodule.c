#include <Python.h>
#include "bgen.h"

typedef struct {
    PyObject_HEAD
    struct bgen_file *bgen;
} BgenReader;

static PyTypeObject BgenReaderType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "simplebgen.BgenReader",             /* tp_name */
    sizeof(BgenReader), /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Noddy objects",           /* tp_doc */
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
    if (PyType_Ready(&BgenReaderType) < 0)
        return NULL;

    m = PyModule_Create(&simplebgenmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&BgenReaderType);
    PyModule_AddObject(m, "BgenReader", (PyObject *)&BgenReaderType);
    return m;
}
