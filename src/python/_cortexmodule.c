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
    dBGraph db_graph;
} Graph;

typedef struct {
    PyObject_HEAD
    Graph *graph;
    uint64_t index;
    char *buffer;
} KmerIterator;

/**********************************************************
 *
 * Graph object 
 *
 **********************************************************/

static void
Graph_dealloc(Graph* self)
{
    db_graph_dealloc(&self->db_graph);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
Graph_init(Graph *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    static char *kwlist[] = {"path", NULL};
    char *path;
    int64_t nkmers;
    GraphFileReader gfile;
    LoadingStats stats = LOAD_STATS_INIT_MACRO;
    GraphLoadingPrefs gprefs = {.db_graph = &self->db_graph,
        .boolean_covgs = false,
        .must_exist_in_graph = false,
        .must_exist_in_edges = NULL,
        .empty_colours = false
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &path)) {
        goto out;
    }
    memset(&gfile, 0, sizeof(gfile));
    if (graph_file_open(&gfile, path) != 1) {
        /* Note that here we have not allocated the db_graph object yet
         * and so we need to put in a flag to detect this. However,
         * all errors result in an abort at the moment, so we can't 
         * do this properly anyway.
         */
        PyErr_SetString(CortexError, "Error opening file");
        goto out;
    }
    nkmers = graph_file_nkmers(&gfile);
    db_graph_alloc(&self->db_graph, gfile.hdr.kmer_size, gfile.hdr.num_of_cols, 
            gfile.hdr.num_of_cols, 2 * nkmers);
    graph_load(&gfile, gprefs, &stats);
    graph_file_close(&gfile);
    ret = 0;
out:
    return ret;
}

static PyObject *
Graph_get_kmer_size(Graph *self)
{
    return Py_BuildValue("I", (unsigned int) self->db_graph.kmer_size); 
}

static PyObject *
Graph_get_num_cols(Graph *self)
{
    return Py_BuildValue("I", (unsigned int) self->db_graph.num_of_cols); 
}

static PyObject *
Graph_get_num_edge_cols(Graph *self)
{
    return Py_BuildValue("I", (unsigned int) self->db_graph.num_edge_cols); 
}

static PyObject *
Graph_get_capacity(Graph *self)
{
    return Py_BuildValue("K", (unsigned long long) self->db_graph.ht.capacity); 
}

static PyMemberDef Graph_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef Graph_methods[] = {
    {"get_kmer_size", (PyCFunction) Graph_get_kmer_size, METH_NOARGS, 
            "Returns the kmer size" },
    {"get_num_cols", (PyCFunction) Graph_get_num_cols, METH_NOARGS, 
            "Returns the number of colours in the graph." },
    {"get_num_edge_cols", (PyCFunction) Graph_get_num_edge_cols, METH_NOARGS, 
            "Returns the number of edge colours in the graph." },
    {"get_capacity", (PyCFunction) Graph_get_capacity, METH_NOARGS, 
            "Returns the capacity of the graph's hash table." },
    {NULL}  /* Sentinel */
};

static PyTypeObject GraphType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_cortex.Graph",             /* tp_name */
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
   
/**********************************************************
 *
 * KmerIterator object 
 *
 **********************************************************/

static void
KmerIterator_dealloc(KmerIterator* self)
{
    PyMem_Free(self->buffer);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int
KmerIterator_init(KmerIterator *self, PyObject *args, PyObject *kwds)
{
    int ret = -1;
    static char *kwlist[] = {"graph", NULL};
    Graph *g; 

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, 
                &GraphType, &g)) {
        goto out;
    }
    self->graph = g;
    self->index = 0;
    self->buffer = PyMem_Malloc(self->graph->db_graph.kmer_size + 1);
    if (self->buffer == NULL) {
        PyErr_NoMemory();
        goto out;
    }
    ret = 0;
out:
    return ret;
}

static PyObject *
KmerIterator_next(KmerIterator *self)
{
    PyObject *ret = NULL;
    BinaryKmer kmer;
    HashTable *hash_table = &self->graph->db_graph.ht;
    size_t k = self->graph->db_graph.kmer_size;
    char *str = self->buffer; 
    while (ret == NULL && self->index < hash_table->capacity) {
        kmer = hash_table->table[self->index];
        if (HASH_ENTRY_ASSIGNED(kmer)) {
            binary_kmer_to_str(kmer, k, str);
            ret = Py_BuildValue("s", str);
            if (ret == NULL) {
                goto out;
            }
        }
        self->index++;
    }
out:
    return ret;
}

static PyMemberDef KmerIterator_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef KmerIterator_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject KmerIteratorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_cortex.KmerIterator",             /* tp_name */
    sizeof(KmerIterator),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)KmerIterator_dealloc, /* tp_dealloc */
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
    "KmerIterator objects",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    PyObject_SelfIter,               /* tp_iter */
    (iternextfunc) KmerIterator_next, /* tp_iternext */
    KmerIterator_methods,             /* tp_methods */
    KmerIterator_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)KmerIterator_init,      /* tp_init */
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
    /* KmerIterator */
    KmerIteratorType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&KmerIteratorType) < 0) {
        INITERROR;
    }
    Py_INCREF(&KmerIteratorType);
    PyModule_AddObject(module, "KmerIterator", (PyObject *) &KmerIteratorType);

    CortexError = PyErr_NewException("_cortex.CortexError", 
            NULL, NULL);
    Py_INCREF(CortexError);
    PyModule_AddObject(module, "CortexError", CortexError);
    

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
