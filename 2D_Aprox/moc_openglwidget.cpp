/****************************************************************************
** Meta object code from reading C++ file 'openglwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "openglwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'openglwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_openglwidget_t {
    QByteArrayData data[17];
    char stringdata0[186];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_openglwidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_openglwidget_t qt_meta_stringdata_openglwidget = {
    {
QT_MOC_LITERAL(0, 0, 12), // "openglwidget"
QT_MOC_LITERAL(1, 13, 9), // "rotate_cw"
QT_MOC_LITERAL(2, 23, 0), // ""
QT_MOC_LITERAL(3, 24, 10), // "rotate_ccw"
QT_MOC_LITERAL(4, 35, 12), // "double_scale"
QT_MOC_LITERAL(5, 48, 9), // "sub_scale"
QT_MOC_LITERAL(6, 58, 9), // "double_nx"
QT_MOC_LITERAL(7, 68, 6), // "sub_nx"
QT_MOC_LITERAL(8, 75, 9), // "double_ny"
QT_MOC_LITERAL(9, 85, 12), // "double_nx_ny"
QT_MOC_LITERAL(10, 98, 9), // "sub_nx_ny"
QT_MOC_LITERAL(11, 108, 6), // "sub_ny"
QT_MOC_LITERAL(12, 115, 11), // "change_func"
QT_MOC_LITERAL(13, 127, 12), // "change_graph"
QT_MOC_LITERAL(14, 140, 13), // "check_threads"
QT_MOC_LITERAL(15, 154, 15), // "increment_error"
QT_MOC_LITERAL(16, 170, 15) // "decrement_error"

    },
    "openglwidget\0rotate_cw\0\0rotate_ccw\0"
    "double_scale\0sub_scale\0double_nx\0"
    "sub_nx\0double_ny\0double_nx_ny\0sub_nx_ny\0"
    "sub_ny\0change_func\0change_graph\0"
    "check_threads\0increment_error\0"
    "decrement_error"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_openglwidget[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   89,    2, 0x0a /* Public */,
       3,    0,   90,    2, 0x0a /* Public */,
       4,    0,   91,    2, 0x0a /* Public */,
       5,    0,   92,    2, 0x0a /* Public */,
       6,    0,   93,    2, 0x0a /* Public */,
       7,    0,   94,    2, 0x0a /* Public */,
       8,    0,   95,    2, 0x0a /* Public */,
       9,    0,   96,    2, 0x0a /* Public */,
      10,    0,   97,    2, 0x0a /* Public */,
      11,    0,   98,    2, 0x0a /* Public */,
      12,    0,   99,    2, 0x0a /* Public */,
      13,    0,  100,    2, 0x0a /* Public */,
      14,    0,  101,    2, 0x0a /* Public */,
      15,    0,  102,    2, 0x0a /* Public */,
      16,    0,  103,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void openglwidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<openglwidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->rotate_cw(); break;
        case 1: _t->rotate_ccw(); break;
        case 2: _t->double_scale(); break;
        case 3: _t->sub_scale(); break;
        case 4: _t->double_nx(); break;
        case 5: _t->sub_nx(); break;
        case 6: _t->double_ny(); break;
        case 7: _t->double_nx_ny(); break;
        case 8: _t->sub_nx_ny(); break;
        case 9: _t->sub_ny(); break;
        case 10: _t->change_func(); break;
        case 11: _t->change_graph(); break;
        case 12: _t->check_threads(); break;
        case 13: _t->increment_error(); break;
        case 14: _t->decrement_error(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject openglwidget::staticMetaObject = { {
    &QOpenGLWidget::staticMetaObject,
    qt_meta_stringdata_openglwidget.data,
    qt_meta_data_openglwidget,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *openglwidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *openglwidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_openglwidget.stringdata0))
        return static_cast<void*>(this);
    return QOpenGLWidget::qt_metacast(_clname);
}

int openglwidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QOpenGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 15)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 15;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
