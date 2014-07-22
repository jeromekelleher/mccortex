#include <Python.h>
#include <structmember.h>

#include "global.h"
#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "assemble_contigs.h"
#include "seq_reader.h"
#include "graph_format.h"
#include "gpath_reader.h"
#include "gpath_checks.h"


#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#define MODULE_DOC \
"Python interface to cortex graphs"

static PyObject *CortexError;

typedef struct {
    PyObject_HEAD
} Graph;


static void
Graph_dealloc(Graph* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}


static int
Graph_init(Graph *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    static char *kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "", kwlist)) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyMemberDef Graph_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef Graph_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject GraphType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_wormtable.Graph",             /* tp_name */
    sizeof(Graph),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Graph_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Graph objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    Graph_methods,             /* tp_methods */
    Graph_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Graph_init,      /* tp_init */
};
   


/* Initialisation code supports Python 2.x and 3.x. The framework uses the 
 * recommended structure from http://docs.python.org/howto/cporting.html. 
 * I've ignored the point about storing state in globals, as the examples 
 * from the Python documentation still use this idiom. 
 */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef cortexmodule = {
    PyModuleDef_HEAD_INIT,
    "_cortex",   /* name of module */
    MODULE_DOC, /* module documentation, may be NULL */
    -1,    
    NULL, NULL, NULL, NULL, NULL 
};

#define INITERROR return NULL

PyObject * 
PyInit__cortex(void)

#else
#define INITERROR return

void
init_cortex(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&cortexmodule);
#else
    PyObject *module = Py_InitModule3("_cortex", NULL, MODULE_DOC);
#endif
    if (module == NULL) {
        INITERROR;
    }
    /* Graph */
    GraphType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&GraphType) < 0) {
        INITERROR;
    }
    Py_INCREF(&GraphType);
    PyModule_AddObject(module, "Graph", (PyObject *) &GraphType);

    CortexError = PyErr_NewException("_cortex.CortexError", 
            NULL, NULL);
    Py_INCREF(CortexError);
    PyModule_AddObject(module, "CortexError", CortexError);
    

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
