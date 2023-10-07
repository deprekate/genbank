#include <Python.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

#if PY_MAJOR >= 3
#define PY3K
#endif

KSEQ_INIT(gzFile, gzread)

typedef struct {
    PyObject_HEAD
	gzFile fp;
	kseq_t *seq;
	int l;
} Sequences;

PyObject* Sequences_iter(PyObject *self){
	Py_INCREF(self);
	return self;
}
PyObject* Sequences_iternext(PyObject *self){
	Sequences *p = (Sequences *)self;
	if((p->l = kseq_read(p->seq)) >= 0){
		return Py_BuildValue("(ss)", 'HEAD' , 'NUCS');
	}else{
		PyErr_SetNone(PyExc_StopIteration);
		return NULL;
	}
}
static void Iter_dealloc(Sequences *self){
	PyObject_Del(self);
	//kseq_destroy(seq);
	//gzclose(fp);
}

static PyTypeObject IterableType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Iter",
    .tp_doc = "Custom objects",
    .tp_basicsize = sizeof(Sequences),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_dealloc = (destructor) Iter_dealloc,
    .tp_iter      = Sequences_iter,
    .tp_iternext  = Sequences_iternext
};
static PyObject* get_sequences(PyObject *self, PyObject *args){ //*Py_UNUSED(ignored)) {
	Sequences *p;
	char *filename;


	p = PyObject_New(Sequences, &IterableType);
	if (!p) return NULL;

	if (!PyObject_Init((PyObject *)p, &IterableType)) {
    	Py_DECREF(p);
    	return NULL;
 	}

	if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

	p->fp = gzopen(filename, "r");
	p->seq = kseq_init(p->fp);

	return 0;
	/* I'm not sure if it's strictly necessary. */
    if (!PyObject_Init((PyObject *)p, &IterableType)) {
        Py_DECREF(p);
        return NULL;
    }
	return (PyObject *) p;
}

// Our Modules Function Definition struct
// We require this `NULL` to signal the end of our method
static PyMethodDef metagenome_methods[] = {
    { "metagenome", (PyCFunction) get_sequences, METH_VARARGS | METH_KEYWORDS, "Returns the edges of connected orfs" },
    { NULL, NULL, 0, NULL }
};
//#ifdef PY3K
// module definition structure for python3
static struct PyModuleDef metagenome = {
     PyModuleDef_HEAD_INIT,
    "metagenome",
    "mod doc",
    -1,
    metagenome_methods
};
PyMODINIT_FUNC PyInit_metagenome(void){
    return PyModule_Create(&metagenome);
}

